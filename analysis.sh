set -eo pipefail
minimap=~/minimap2-2.16_x64-linux/minimap2
hs37d5=/research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa
sentieon_hapA=/research/rv-02/home/eanderson/Sentieon_Stratification_Chr19/synthetic_binned_fasta/sentieon_pe400_chr19/hapA/binned.fasta
sentieon_hapB=/research/rv-02/home/eanderson/Sentieon_Stratification_Chr19/synthetic_binned_fasta/sentieon_pe400_chr19/hapB/binned.fasta
in_house_hapA=/home/ycai/hm/dev/in_house_denovo/phasing/new/chr19/seed_more/snp7_20k/merge_hom/spades_hapa/contigs.fasta
in_house_hapB=/home/ycai/hm/dev/in_house_denovo/phasing/new/chr19/seed_more/snp7_20k/merge_hom/spades_hapb/contigs.fasta
truth_hapA=/research/rv-02/home/eanderson/Simulated_Data_Wenlan/Simulated_Ref/Fasta/chr19_normalized_hapA.fa
truth_hapB=/research/rv-02/home/eanderson/Simulated_Data_Wenlan/Simulated_Ref/Fasta/chr19_normalized_hapB.fa

mkdir -p unique_aligns_to_hs37d5
$minimap -x asm5 -t 20 --secondary=no --cs -y $hs37d5 <(sed 's/=/:i:/g' $sentieon_hapA) > unique_aligns_to_hs37d5/sentieon_hapA.paf
$minimap -x asm5 -t 20 --secondary=no --cs -y $hs37d5 <(sed 's/=/:i:/g' $sentieon_hapB) > unique_aligns_to_hs37d5/sentieon_hapB.paf
$minimap -x asm5 -t 20 --secondary=no --cs -y $hs37d5 <(sed 's/=/:i:/g' $in_house_hapA) > unique_aligns_to_hs37d5/in_house_hapA.paf
$minimap -x asm5 -t 20 --secondary=no --cs -y $hs37d5 <(sed 's/=/:i:/g' $in_house_hapB) > unique_aligns_to_hs37d5/in_house_hapB.paf
$minimap -x asm5 -t 20 --secondary=no --cs -y $hs37d5 <(sed 's/=/:i:/g' $truth_hapA) > unique_aligns_to_hs37d5/truth_hapA.paf
$minimap -x asm5 -t 20 --secondary=no --cs -y $hs37d5 <(sed 's/=/:i:/g' $truth_hapB) > unique_aligns_to_hs37d5/truth_hapB.paf

paftools.js call <(cat unique_aligns_to_hs37d5/sentieon_hapA.paf | sort -k6,6 -k8,8n) > unique_aligns_to_hs37d5/sentieon_hapA.var
paftools.js call <(cat unique_aligns_to_hs37d5/sentieon_hapB.paf | sort -k6,6 -k8,8n) > unique_aligns_to_hs37d5/sentieon_hapB.var
paftools.js call <(cat unique_aligns_to_hs37d5/in_house_hapA.paf | sort -k6,6 -k8,8n) > unique_aligns_to_hs37d5/in_house_hapA.var
paftools.js call <(cat unique_aligns_to_hs37d5/in_house_hapB.paf | sort -k6,6 -k8,8n) > unique_aligns_to_hs37d5/in_house_hapB.var
paftools.js call <(cat unique_aligns_to_hs37d5/truth_hapA.paf | sort -k6,6 -k8,8n) > unique_aligns_to_hs37d5/truth_hapA.var
paftools.js call <(cat unique_aligns_to_hs37d5/truth_hapB.paf | sort -k6,6 -k8,8n) > unique_aligns_to_hs37d5/truth_hapB.var



grep ^R unique_aligns_to_hs37d5/sentieon_hapA.var | awk '{sum+=$4-$3} END {print (sum)}' > unique_aligns_to_hs37d5/sentieon_hapA_unique_align_length.txt
grep ^R unique_aligns_to_hs37d5/sentieon_hapB.var | awk '{sum+=$4-$3} END {print (sum)}' > unique_aligns_to_hs37d5/sentieon_hapB_unique_align_length.txt
grep ^R unique_aligns_to_hs37d5/in_house_hapA.var | awk '{sum+=$4-$3} END {print (sum)}' > unique_aligns_to_hs37d5/in_house_hapA_unique_align_length.txt
grep ^R unique_aligns_to_hs37d5/in_house_hapB.var | awk '{sum+=$4-$3} END {print (sum)}' > unique_aligns_to_hs37d5/in_house_hapB_unique_align_length.txt
grep ^R unique_aligns_to_hs37d5/truth_hapA.var | awk '{sum+=$4-$3} END {print (sum)}' > unique_aligns_to_hs37d5/truth_hapA_unique_align_length.txt
grep ^R unique_aligns_to_hs37d5/truth_hapB.var | awk '{sum+=$4-$3} END {print (sum)}' > unique_aligns_to_hs37d5/truth_hapB_unique_align_length.txt

grep ^R unique_aligns_to_hs37d5/sentieon_hapA.var | awk '{print $2"\t"$3"\t"$4}' > unique_aligns_to_hs37d5/sentieon_hapA_unique_aligns.bed
grep ^R unique_aligns_to_hs37d5/sentieon_hapB.var | awk '{print $2"\t"$3"\t"$4}' > unique_aligns_to_hs37d5/sentieon_hapB_unique_aligns.bed
grep ^R unique_aligns_to_hs37d5/in_house_hapA.var | awk '{print $2"\t"$3"\t"$4}' > unique_aligns_to_hs37d5/in_house_hapA_unique_aligns.bed
grep ^R unique_aligns_to_hs37d5/in_house_hapB.var | awk '{print $2"\t"$3"\t"$4}' > unique_aligns_to_hs37d5/in_house_hapB_unique_aligns.bed
grep ^R unique_aligns_to_hs37d5/truth_hapA.var | awk '{print $2"\t"$3"\t"$4}' > unique_aligns_to_hs37d5/truth_hapA_unique_aligns.bed
grep ^R unique_aligns_to_hs37d5/truth_hapB.var | awk '{print $2"\t"$3"\t"$4}' > unique_aligns_to_hs37d5/truth_hapB_unique_aligns.bed

mkdir -p all_alignments_to_hs37d5
bedtools sort -i <(cat unique_aligns_to_hs37d5/sentieon_hapA.paf | cut -f 6,8,9) -g ${hs37d5}.fai > all_alignments_to_hs37d5/sentieon_hapA.bed
bedtools sort -i <(cat unique_aligns_to_hs37d5/sentieon_hapB.paf | cut -f 6,8,9) -g ${hs37d5}.fai > all_alignments_to_hs37d5/sentieon_hapB.bed
bedtools sort -i <(cat unique_aligns_to_hs37d5/in_house_hapA.paf | cut -f 6,8,9) -g ${hs37d5}.fai > all_alignments_to_hs37d5/in_house_hapA.bed
bedtools sort -i <(cat unique_aligns_to_hs37d5/in_house_hapB.paf | cut -f 6,8,9) -g ${hs37d5}.fai > all_alignments_to_hs37d5/in_house_hapB.bed
bedtools sort -i <(cat unique_aligns_to_hs37d5/truth_hapA.paf | cut -f 6,8,9) -g ${hs37d5}.fai > all_alignments_to_hs37d5/truth_hapA.bed
bedtools sort -i <(cat unique_aligns_to_hs37d5/truth_hapB.paf | cut -f 6,8,9) -g ${hs37d5}.fai > all_alignments_to_hs37d5/truth_hapB.bed

awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/sentieon_hapA.bed > all_alignments_to_hs37d5/sentieon_hapA_total_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/sentieon_hapB.bed > all_alignments_to_hs37d5/sentieon_hapB_total_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/in_house_hapA.bed > all_alignments_to_hs37d5/in_house_hapA_total_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/in_house_hapB.bed > all_alignments_to_hs37d5/in_house_hapB_total_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/truth_hapA.bed > all_alignments_to_hs37d5/truth_hapA_total_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/truth_hapB.bed > all_alignments_to_hs37d5/truth_hapB_total_coverage.txt

bedtools merge -i all_alignments_to_hs37d5/sentieon_hapA.bed > all_alignments_to_hs37d5/sentieon_hapA_compressed_coverage.bed
bedtools merge -i all_alignments_to_hs37d5/sentieon_hapB.bed > all_alignments_to_hs37d5/sentieon_hapB_compressed_coverage.bed
bedtools merge -i all_alignments_to_hs37d5/in_house_hapA.bed > all_alignments_to_hs37d5/in_house_hapA_compressed_coverage.bed
bedtools merge -i all_alignments_to_hs37d5/in_house_hapB.bed > all_alignments_to_hs37d5/in_house_hapB_compressed_coverage.bed
bedtools merge -i all_alignments_to_hs37d5/truth_hapA.bed > all_alignments_to_hs37d5/truth_hapA_compressed_coverage.bed
bedtools merge -i all_alignments_to_hs37d5/truth_hapB.bed > all_alignments_to_hs37d5/truth_hapB_compressed_coverage.bed

awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/sentieon_hapA_compressed_coverage.bed > all_alignments_to_hs37d5/sentieon_hapA_compressed_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/sentieon_hapB_compressed_coverage.bed > all_alignments_to_hs37d5/sentieon_hapB_compressed_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/in_house_hapA_compressed_coverage.bed > all_alignments_to_hs37d5/in_house_hapA_compressed_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/in_house_hapB_compressed_coverage.bed > all_alignments_to_hs37d5/in_house_hapB_compressed_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/truth_hapA_compressed_coverage.bed > all_alignments_to_hs37d5/truth_hapA_compressed_coverage.txt
awk '{sum+=$3-$2} END {print sum}' all_alignments_to_hs37d5/truth_hapB_compressed_coverage.bed > all_alignments_to_hs37d5/truth_hapB_compressed_coverage.txt

mkdir -p multi_aligns_to_hs37d5
bedtools subtract -a all_alignments_to_hs37d5/sentieon_hapA_compressed_coverage.bed -b unique_aligns_to_hs37d5/sentieon_hapA_unique_aligns.bed > multi_aligns_to_hs37d5/sentieon_hapA_multi_aligns.bed
bedtools subtract -a all_alignments_to_hs37d5/sentieon_hapB_compressed_coverage.bed -b unique_aligns_to_hs37d5/sentieon_hapB_unique_aligns.bed > multi_aligns_to_hs37d5/sentieon_hapB_multi_aligns.bed
bedtools subtract -a all_alignments_to_hs37d5/in_house_hapA_compressed_coverage.bed -b unique_aligns_to_hs37d5/in_house_hapA_unique_aligns.bed > multi_aligns_to_hs37d5/in_house_hapA_multi_aligns.bed
bedtools subtract -a all_alignments_to_hs37d5/in_house_hapB_compressed_coverage.bed -b unique_aligns_to_hs37d5/in_house_hapB_unique_aligns.bed > multi_aligns_to_hs37d5/in_house_hapB_multi_aligns.bed
bedtools subtract -a all_alignments_to_hs37d5/truth_hapA_compressed_coverage.bed -b unique_aligns_to_hs37d5/truth_hapA_unique_aligns.bed > multi_aligns_to_hs37d5/truth_hapA_multi_aligns.bed
bedtools subtract -a all_alignments_to_hs37d5/truth_hapB_compressed_coverage.bed -b unique_aligns_to_hs37d5/truth_hapB_unique_aligns.bed > multi_aligns_to_hs37d5/truth_hapB_multi_aligns.bed

awk '{sum+=$3-$2} END {print sum}' multi_aligns_to_hs37d5/sentieon_hapA_multi_aligns.bed > multi_aligns_to_hs37d5/sentieon_hapA_multi_aligns_coverage.txt
awk '{sum+=$3-$2} END {print sum}' multi_aligns_to_hs37d5/sentieon_hapB_multi_aligns.bed > multi_aligns_to_hs37d5/sentieon_hapB_multi_aligns_coverage.txt
awk '{sum+=$3-$2} END {print sum}' multi_aligns_to_hs37d5/in_house_hapA_multi_aligns.bed > multi_aligns_to_hs37d5/in_house_hapA_multi_aligns_coverage.txt
awk '{sum+=$3-$2} END {print sum}' multi_aligns_to_hs37d5/in_house_hapB_multi_aligns.bed > multi_aligns_to_hs37d5/in_house_hapB_multi_aligns_coverage.txt
awk '{sum+=$3-$2} END {print sum}' multi_aligns_to_hs37d5/truth_hapA_multi_aligns.bed > multi_aligns_to_hs37d5/truth_hapA_multi_aligns_coverage.txt
awk '{sum+=$3-$2} END {print sum}' multi_aligns_to_hs37d5/truth_hapB_multi_aligns.bed > multi_aligns_to_hs37d5/truth_hapB_multi_aligns_coverage.txt

bedtools genomecov -bg -i all_alignments_to_hs37d5/sentieon_hapA.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai > all_alignments_to_hs37d5/sentieon_hapA.genome_cov
bedtools genomecov -bg -i all_alignments_to_hs37d5/sentieon_hapB.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai > all_alignments_to_hs37d5/sentieon_hapB.genome_cov
bedtools genomecov -bg -i all_alignments_to_hs37d5/in_house_hapA.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai > all_alignments_to_hs37d5/in_house_hapA.genome_cov
bedtools genomecov -bg -i all_alignments_to_hs37d5/in_house_hapB.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai > all_alignments_to_hs37d5/in_house_hapB.genome_cov
bedtools genomecov -bg -i all_alignments_to_hs37d5/truth_hapA.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai > all_alignments_to_hs37d5/truth_hapA.genome_cov
bedtools genomecov -bg -i all_alignments_to_hs37d5/truth_hapB.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai > all_alignments_to_hs37d5/truth_hapB.genome_cov

bedtools genomecov -d -i all_alignments_to_hs37d5/sentieon_hapA.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai | awk '{if( $3 > 0){print $0}}' > multi_aligns_to_hs37d5/sentieon_hapA.genome_cov
bedtools genomecov -d -i all_alignments_to_hs37d5/sentieon_hapB.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai | awk '{if( $3 > 0){print $0}}' > multi_aligns_to_hs37d5/sentieon_hapB.genome_cov
bedtools genomecov -d -i all_alignments_to_hs37d5/in_house_hapA.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai | awk '{if( $3 > 0){print $0}}' > multi_aligns_to_hs37d5/in_house_hapA.genome_cov
bedtools genomecov -d -i all_alignments_to_hs37d5/in_house_hapB.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai | awk '{if( $3 > 0){print $0}}' > multi_aligns_to_hs37d5/in_house_hapB.genome_cov
bedtools genomecov -d -i all_alignments_to_hs37d5/truth_hapA.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai | awk '{if( $3 > 0){print $0}}' > multi_aligns_to_hs37d5/truth_hapA.genome_cov
bedtools genomecov -d -i all_alignments_to_hs37d5/truth_hapB.bed -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai | awk '{if( $3 > 0){print $0}}' > multi_aligns_to_hs37d5/truth_hapB.genome_cov


mkdir -p Unaligned_Truth_Regions
bedtools intersect -a all_alignments_to_hs37d5/truth_hapA_compressed_coverage.bed -b all_alignments_to_hs37d5/truth_hapB_compressed_coverage.bed > Unaligned_Truth_Regions/diploid_truth_regions.bed

bedtools merge -i \
      <(bedtools sort -i <(cat all_alignments_to_hs37d5/in_house_hapA_compressed_coverage.bed \
                           all_alignments_to_hs37d5/in_house_hapB_compressed_coverage.bed \
                           all_alignments_to_hs37d5/sentieon_hapA_compressed_coverage.bed \
                           all_alignments_to_hs37d5/sentieon_hapB_compressed_coverage.bed) \
                      -g ${hs37d5}.fai) > Unaligned_Truth_Regions/all_in_house_and_sentieon_coverage.bed

bedtools merge -i \
      <(bedtools sort -i <(cat all_alignments_to_hs37d5/in_house_hapA_compressed_coverage.bed \
                           all_alignments_to_hs37d5/in_house_hapB_compressed_coverage.bed) \
                      -g ${hs37d5}.fai) > Unaligned_Truth_Regions/all_in_house_coverage.bed

bedtools merge -i \
    <(bedtools sort -i <(cat all_alignments_to_hs37d5/sentieon_hapA_compressed_coverage.bed \
                         all_alignments_to_hs37d5/sentieon_hapB_compressed_coverage.bed) \
                    -g ${hs37d5}.fai) > Unaligned_Truth_Regions/all_sentieon_coverage.bed

bedtools intersect -a Unaligned_Truth_Regions/all_sentieon_coverage.bed \
                   -b Unaligned_Truth_Regions/all_in_house_coverage.bed > Unaligned_Truth_Regions/all_shared_coverage.bed

bedtools subtract -a Unaligned_Truth_Regions/diploid_truth_regions.bed \
                  -b Unaligned_Truth_Regions/all_in_house_and_sentieon_coverage.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns.bed

bedtools subtract -a Unaligned_Truth_Regions/diploid_truth_regions.bed \
                  -b Unaligned_Truth_Regions/all_sentieon_coverage.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed

bedtools subtract -a Unaligned_Truth_Regions/diploid_truth_regions.bed \
                  -b Unaligned_Truth_Regions/all_in_house_coverage.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed

bedtools subtract -a Unaligned_Truth_Regions/diploid_truth_regions.bed \
                  -b Unaligned_Truth_Regions/all_shared_coverage.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_shared.bed

bedtools slop -b 1000 -i Unaligned_Truth_Regions/diploid_truth_regions_no_aligns.bed -g ${hs37d5}.fai > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_1kb_slop.bed 
bedtools slop -b 1000 -i Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed -g ${hs37d5}.fai > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon_1kb_slop.bed

awk '{sum+=$3-$2} END {print sum, sum/NR}' Unaligned_Truth_Regions/diploid_truth_regions_no_aligns.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns.txt
awk '{sum+=$3-$2} END {print sum, sum/NR}' Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.txt
awk '{sum+=$3-$2} END {print sum, sum/NR}' Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.txt

python3 ../get_variants_for_bed_test.py Unaligned_Truth_Regions/diploid_truth_regions_no_aligns.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_het_vars.bed
python3 ../get_variants_for_bed_test.py Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house_het_vars.bed
python3 ../get_variants_for_bed_test.py Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon_het_vars.bed

python3 ../get_variants_for_bed_test.py Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_1kb_slop.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_het_vars_1kb_slop.bed
python3 ../get_variants_for_bed_test.py Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon_1kb_slop.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon_het_vars_1kb_slop.bed

awk '{if($6 > 0.5) { print $0}}' Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house_het_vars.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house_high_het.bed
awk '{if($6 <= 0.5) { print $0}}' Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house_het_vars.bed > Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house_low_het.bed

