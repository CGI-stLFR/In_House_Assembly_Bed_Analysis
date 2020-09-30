import os

configfile: "config.yaml"

# function to define target files
def run_all_files():
    expected_files = []
    for ref in config['ref_files']:
        for sample in config['samples']:
            outfile = "Difficult_Regions_Summary/" + '/' + ref + '/' + sample + "/summary.txt"
            expected_files.append(outfile)
            for hap in ['0', '1']:
                outfile = "Genome_Coverage_Summary/" + ref + '/' + sample + '/' + hap + "/genome_coverage_summary.txt"
                expected_files.append(outfile)
                outfile = "Aligns_as_Bed/" + ref + '/' + sample + '/' + hap + "/ref_aligns.bed"
                expected_files.append(outfile)

            for region in config['difficult_regions']:
                outfile = "Stratified_Beds_No_Aligns/" + ref + '/' + sample + '/' + region + "/missing_difficult_regions.bed"
                expected_files.append(outfile)
                outfile = "Stratified_Beds_With_Aligns/" + ref + '/' + sample + '/' + region + "/covered_difficult_regions.bed"
                expected_files.append(outfile)


            outfile = "Heterozygosity_Summary/" + ref  + '/' + sample + "/heterozygosity_summary.txt"
            expected_files.append(outfile)

        for truth_hap in ['hapA', 'hapB']:
            outfile = "Truth_Genome_Coverage_Summary/" + ref + '/' + truth_hap + "/genome_coverage_summary.txt"
            expected_files.append(outfile)

        for region in config['difficult_regions']:
            outfile = "Diploid_Difficult_Regions/" + ref + '/' + region + "diploid_difficult_regions.bed"


# rule to generate targets
rule run_all:
    input:
        all_outputs = run_all_files()


# run minimap on all samples and haplotypes
rule run_minimap_on_samples:
    input:
        fasta = "Assembly/{sample}/contigs_{hap}.fa",
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    output:
        "Align_to_Ref/{ref}/{sample}/{hap}/ref_aligns.paf"
    params:
        minimap = config['minimap2']
    threads:
        config['threads']['minimap']
    shell:
        "{params.minimap} -x asm5 "
            "-t {threads} "
            "--secondary=no "
            "--cs "
            "-y "
            "{input.ref} "
            "<( sed 's/=/:i:/g' {input.fasta}) "
            "> {output}"


# run minimap on truth samples
rule run_minimap_on_truths:
    input:
        fasta = config['truth_dir'] + "chr19_normalized_{truth_hap}.fa",
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    output:
        "Align_truth_to_Ref/{ref}/{truth_hap}/ref_aligns.paf"
    params:
        minimap = config['minimap2']
    threads:
        config['threads']['minimap']
    shell:
        "{params.minimap} -x asm5 "
            "-t {threads} "
            "--secondary=no "
            "--cs "
            "-y "
            "{input.ref} "
            "<( sed 's/=/:i:/g' {input.fasta}) "
            "> {output}"


# get bed file from paf files
rule get_bed_file_from_paf:
    input:
        paf = "Align_to_Ref/{ref}/{sample}/{hap}/ref_aligns.paf",
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    output:
        "Aligns_as_Bed/{ref}/{sample}/{hap}/ref_aligns.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} sort -i <(cat {input.paf} | cut -f 6,8,9) -g {input.ref}.fai > {output}"


# get bed file from truth pafs
rule get_truth_bed_file_from_paf:
    input:
        paf = "Align_truth_to_Ref/{ref}/{truth_hap}/ref_aligns.paf",
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    output:
        "Truth_Aligns_as_Bed/{ref}/{truth_hap}/ref_aligns.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} sort -i <(cat {input.paf} | cut -f 6,8,9) -g {input.ref}.fai > {output}"


# get a genomecoverage bedgraph from paf beds
rule create_bedgraph_genomecov:
    input:
        bed = "Aligns_as_Bed/{ref}/{sample}/{hap}/ref_aligns.bed",
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    output:
        "Genome_Coverage/{ref}/{sample}/{hap}/ref_aligns.genome_cov"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} genomecov -bg -i {input.bed} -g {input.ref}.fai > {output}"


# get a genomecoverage bedgraph from truth paf beds
rule create_truth_bedgraph_genomecov:
    input:
        bed = "Truth_Aligns_as_Bed/{ref}/{truth_hap}/ref_aligns.bed",
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    output:
        "Truth_Genome_Coverage/{ref}/{truth_hap}/ref_aligns.genome_cov"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} genomecov -bg -i {input.bed} -g {input.ref}.fai > {output}"


# process genomecoverage beds to get various coverages and summarize
rule process_genomecov_files:
    input:
        "Genome_Coverage/{ref}/{sample}/{hap}/ref_aligns.genome_cov"
    output:
        "Genome_Coverage_Summary/{ref}/{sample}/{hap}/genome_coverage_summary.txt"
    run:
        unique_cov = 0
        multi_cov = 0
        with open(input[0], 'r') as gc, open(output[0], 'w') as outf:
            for line in gc:
                parts = line.strip().split()
                if int(parts[3]) == 1:
                    unique_cov += int(parts[2]) - int(parts[1]) - 1
                else:
                    multi_cov += (int(parts[2]) - int(parts[1]) - (1 * int(parts[3]))) * int(parts[3])

            print(f"Unique Coverage:\t{unique_cov}\nMulti Coverage:\t{multi_cov}\nTotal Coverage:\t{multi_cov + unique_cov}", file=outf)


# process truth genomecoverage beds to get summary
rule process_truth_genomecov_files:
    input:
        "Truth_Genome_Coverage/{ref}/{truth_hap}/ref_aligns.genome_cov"
    output:
        "Truth_Genome_Coverage_Summary/{ref}/{truth_hap}/genome_coverage_summary.txt"
    run:
        unique_cov = 0
        multi_cov = 0
        with open(input[0], 'r') as gc, open(output[0], 'w') as outf:
            for line in gc:
                parts = line.strip().split()
                if int(parts[3]) == 1:
                    unique_cov += int(parts[2]) - int(parts[1]) - 1
                else:
                    multi_cov += (int(parts[2]) - int(parts[1]) - (1 * int(parts[3]))) * int(parts[3])

            print(f"Unique Coverage:\t{unique_cov}\nMulti Coverage:\t{multi_cov}\nTotal Coverage:\t{multi_cov + unique_cov}", file=outf)


# get the diploid coverage area of the truth haplotypes
rule get_truth_diploid_coverage:
    input:
        expand("Truth_Aligns_as_Bed/{ref}/{truth_hap}/ref_aligns.bed", ref="{ref}", truth_hap=["hapA", "hapB"])
    output:
        "Truth_Diploid_Coverage/{ref}/diploid_truth_regions.bed"
    params:
        bedtools = config['bedtools']
    run:
        command = [params.bedtools, "intersect", "-a", input[0], "-b", input[1],
                   ">", output[0]]
        shell(" ".join(command))


# get combined coverage of sample beds in one file
rule get_combined_coverage:
    input:
        beds = expand("Aligns_as_Bed/{ref}/{sample}/{hap}/ref_aligns.bed", ref="{ref}", sample="{sample}", hap=['0','1']),
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    output:
        "Combined_Coverage/{ref}/{sample}/combined_coverage.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} merge -i "
            "<({params.bedtools} sort -i "
                "<(cat {input.beds} ) "
                "-g {input.ref}.fai ) "
            "> {output}"


# get diploid truth regions without alignments of various samples
rule get_truth_regions_without_aligns:
    input:
        bed = "Combined_Coverage/{ref}/{sample}/combined_coverage.bed",
        truth = "Truth_Diploid_Coverage/{ref}/diploid_truth_regions.bed"
    output:
        "Unaligned_Truth_Regions/{ref}/{sample}/diploid_truth_regions_no_aligns.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} subtract -a {input.truth} -b {input.bed} > {output}"


# annotate the truth regions without alignments with variant data
rule get_variant_annotated_bed:
    input:
        "Unaligned_Truth_Regions/{ref}/{sample}/diploid_truth_regions_no_aligns.bed"
    output:
        "Unaligned_Truth_Vars/{ref}/{sample}/diploid_truth_regions_no_aligns_w_variants.bed"
    params:
        get_variants = config['variants_script'],
        truth_vcf = config['truth_vcf']
    shell:
        "python3 {params.get_variants} {input} {params.truth_vcf} > {output}"


# parse the truth regions annotations and summarize
rule parse_variant_bed:
    input:
        "Unaligned_Truth_Vars/{ref}/{sample}/diploid_truth_regions_no_aligns_w_variants.bed"
    output:
        "Heterozygosity_Summary/{ref}/{sample}/heterozygosity_summary.txt"
    run:
        import pandas as pd

        with open(output[0], 'w') as outf:

            df = pd.read_csv(input[0], sep='\t', names = ['Chr', 'Start', 'End', 'Len', 'NVars', 'VarsPerKB'])
            print("mean length:\t", df['Len'].mean(), file=outf)
            print("total missing:\t", df['Len'].sum(), file=outf)

            tmp_bed_regions = df[df['VarsPerKB'] > 0.5]
            mean=tmp_bed_regions['Len'].mean()
            print("Sum absent Coverage > 0.5:\t", tmp_bed_regions['Len'].sum(), file=outf)
            print("Percent length:\t", tmp_bed_regions['Len'].sum()/df['Len'].sum(), file=outf)
            print("Mean Length:\t", mean, file=outf)

            tmp_bed_regions = df[df['VarsPerKB'] <= 0.5]
            mean=tmp_bed_regions['Len'].mean()
            print("Sum absent Coverage <= 0.5:\t", tmp_bed_regions['Len'].sum(), file=outf)
            print("Percent length:\t", tmp_bed_regions['Len'].sum()/df['Len'].sum(), file=outf)
            print("Mean Length:\t", mean, file=outf)


# return difficult region based on wildcard
def get_difficult_regions_bed(wildcards):
    return config['difficult_regions'][wildcards.region]


# stratify the unaligned regions
rule stratify_difficult_regions:
    input:
        "Unaligned_Truth_Regions/{ref}/{sample}/diploid_truth_regions_no_aligns.bed"
    output:
        "Stratified_Beds_No_Aligns/{ref}/{sample}/{region}/missing_difficult_regions.bed"
    params:
        region = get_difficult_regions_bed,
        bedtools = config['bedtools']
    shell:
        "bedtools intersect -a {input} -b {params.region} > {output}"


# stratify truth regions
rule stratify_truth_with_difficult_regions:
    input:
        "Truth_Diploid_Coverage/{ref}/diploid_truth_regions.bed"
    output:
        "Diploid_Difficult_Regions/{ref}/{region}/diploid_difficult_regions.bed"
    params:
        region = get_difficult_regions_bed,
        bedtools = config['bedtools']
    shell:
        "bedtools intersect -a {input} -b {params.region} > {output}"


# stratify difficult regions with coverage
rule get_covered_difficult_regions:
    input:
        truth = "Diploid_Difficult_Regions/{ref}/{region}/diploid_difficult_regions.bed",
        bed = "Combined_Coverage/{ref}/{sample}/combined_coverage.bed"
    output:
        "Stratified_Beds_With_Aligns/{ref}/{sample}/{region}/covered_difficult_regions.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} intersect -a {input.bed} -b {input.truth} > {output}"


# combine all difficult regions without coverage
rule get_total_difficult_regions_no_cov:
    input:
        expand("Stratified_Beds/{ref}/{sample}/{region}/difficult_regions_no_aligns.bed", ref="{ref}", sample="{sample}", regions=config['difficult_regions'].keys())
    output:
        "Combined_Stratified_Regions_No_Aligns/{ref}/{sample}/combined_no_aligns.bed"
    params:
        bedtools = config['bedtools'],
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    shell:
        "{params.bedtools} merge -i "
            "<({params.bedtools} sort -i <( cat {input} ) "
            "-g {params.ref}.fai) > {output}"


# combine all difficult regions with coverage
rule get_total_difficult_regions_w_cov:
    input:
        expand("Covered_Difficult_Truth_Regions/{ref}/{sample}/{region}/covered_difficult_regions.bed", ref="{ref}", sample="{sample}", regions=config['difficult_regions'].keys())
    output:
        "Combined_Covered_Stratified_Regions/{ref}/{sample}/combined_covered.bed"
    params:
        bedtools = config['bedtools'],
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    shell:
        "{params.bedtools} merge -i "
            "<({params.bedtools} sort -i <( cat {input} ) "
            "-g {params.ref}.fai) > {output}"


# parse individual and aggregate difficult regions with and without coverage
# create a summary file
rule parse_difficult_regions:
    input:
        ind_covered = expand("Stratified_Beds_With_Aligns/{ref}/{sample}/{region}/covered_difficult_regions.bed", ref="{ref}", sample="{sample}", regions=config['difficult_regions'].keys()),
        ind_missing = expand("Stratified_Beds_No_Aligns/{ref}/{sample}/{region}/missing_difficult_regions.bed", ref="{ref}", sample="{sample}", regions=config['difficult_regions'].keys()),
        tot_covered = "Combined_Covered_Stratified_Regions/{ref}/{sample}/combined_covered.bed",
        tot_missing = "Combined_Stratified_Regions_No_Aligns/{ref}/{sample}/combined_no_aligns.bed"
    output:
        "Difficult_Regions_Summary/{ref}/{sample}/summary.txt"
    run:
        regions = config['difficult_regions'].keys()


        ind_covered_results = {}
        for bedfile in ind_covered:
            bed_total = 0
            region = bedfile.split("/")[3]
            with open(bedfile, "r") as bed:
                for line in bed:
                    parts = line.split()
                    bed_total = int(parts[2]) - int(parts[1]) - 1

            ind_covered_results[region] = bed_total


        ind_missing_results = {}
        for bedfile in ind_missing:
            bed_total = 0
            region = bedfile.split("/")[3]
            with open(bedfile, "r") as bed:
                for line in bed:
                    parts = line.split()
                    bed_total = int(parts[2]) - int(parts[1]) - 1

            ind_missing_results[region] = bed_total


        tot_covered_results = 0
        with open(tot_covered, "r") as bed:
            for line in bed:
                parts = line.split()
                tot_covered_results = int(parts[2]) - int(parts[1]) - 1


        tot_missing_results = 0
        with open(tot_missing, "r") as bed:
            for line in bed:
                parts = line.split()
                tot_missing_results = int(parts[2]) - int(parts[1]) - 1

        with open(outfile, "w") as outf:
            print("Total Covered Difficult Regions:", file=outf)
            print(tot_covered_results, file=outf)

            print("Individual Covered Regions:", file=outf)
            for region in regions:
                print(region + ":", file=outf)
                print(ind_covered_results[region], file=outf)

            print("Total Missing Difficult Regions:", file=outf)
            print(tot_missing_results, file=outf)

            print("Individual Missing Regions:", file=outf)
            for region in regions:
                print(region + ":", file=outf)
                print(ind_missing_results[region], file=outf)
