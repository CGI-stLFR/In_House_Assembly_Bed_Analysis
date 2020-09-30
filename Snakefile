import os

configfile: "config.yaml"

def run_all_files():
    expected_files = []
    for ref in config['ref_files']:
        for sample in config['samples']:
            for hap in ['0', '1']:
                file = "Genome_Coverage_Summary/" + ref + '/' + sample + '/' + hap + "/genome_coverage_summary.txt"
                expected_files.append(file)
                file = "Aligns_as_Bed/" + ref + '/' + sample + '/' + hap + "/ref_aligns.bed"
                expected_files.append(file)

            file = "Heterozygosity_Summary/" + ref  + '/' + sample + "/heterozygosity_summary.txt"
            expected_files.append(file)

        for truth_hap in ['hapA', 'hapB']:
            file = "Truth_Genome_Coverage_Summary/" + ref + '/' + truth_hap + "/genome_coverage_summary.txt"
            expected_files.append(file)

    return expected_files


rule run_all:
    input:
        all_outputs = run_all_files()


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


rule compress_coverage:
    input:
        "Aligns_as_Bed/{ref}/{sample}/{hap}/ref_aligns.bed"
    output:
        "Compressed_Coverage/{ref}/{sample}/{hap}/compressed_coverage.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} merge -i {input} > {output}"


rule compress_truth_coverage:
    input:
        "Truth_Aligns_as_Bed/{ref}/{truth_hap}/ref_aligns.bed"
    output:
        "Truth_Compressed_Coverage/{ref}/{truth_hap}/compressed_coverage.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} merge -i {input} > {output}"


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

rule get_combined_coverage:
    input:
        beds = expand("Aligns_as_Bed/{ref}/{sample}/{hap}/ref_aligns.bed", ref="{ref}", sample="{sample}", hap=['0','1']),
        ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref])
    output:
        "Diploid_Coverage/{ref}/{sample}/diploid_coverage.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} merge -i "
            "<({params.bedtools} sort -i "
                "<(cat {input.beds} ) "
                "-g {input.ref}.fai ) "
            "> {output}"


rule get_truth_regions_without_aligns:
    input:
        bed = "Diploid_Coverage/{ref}/{sample}/diploid_coverage.bed",
        truth = "Truth_Diploid_Coverage/{ref}/diploid_truth_regions.bed"
    output:
        "Unaligned_Truth_Regions/{ref}/{sample}/diploid_truth_regions_no_aligns.bed"
    params:
        bedtools = config['bedtools']
    shell:
        "{params.bedtools} subtract -a {input.truth} -b {input.bed} > {output}"


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


