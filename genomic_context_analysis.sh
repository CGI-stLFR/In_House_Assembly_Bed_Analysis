# get intersection of bed files
bedtools intersect -a Bed_Files/LowComplexity/Human_Full_Genome_TRDB_hg19_150331_all_merged_coordonly.bed \
    -b ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed > \
    no_aligns_in_house_lowcmp.bed
bedtools intersect -a Bed_Files/LowComplexity/Human_Full_Genome_TRDB_hg19_150331_all_merged_coordonly.bed \
    -b ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed > \
    no_aligns_sentieon_lowcmp.bed

bedtools intersect -a Bed_Files/SegmentalDuplications/hg19_self_chain_split_both_coordonly.bed \
    -b ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed > \
    no_aligns_in_house_segdup.bed
bedtools intersect -a Bed_Files/SegmentalDuplications/hg19_self_chain_split_both_coordonly.bed \
    -b ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed > \
    no_aligns_sentieon_segdup.bed

bedtools intersect -a Bed_Files/mappability/human_g1k_v37_gemmap_l100_m0_e0_nonuniq.sort_coordonly.bed \
    -b ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed > \
    no_aligns_in_house_lowmap.bed
bedtools intersect -a Bed_Files/mappability/human_g1k_v37_gemmap_l100_m0_e0_nonuniq.sort_coordonly.bed \
    -b ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed > \
    no_aligns_sentieon_lowmap.bed

# get sum of regions
awk '{sum += $3-$2-1} END {print sum}' no_aligns_in_house_lowcmp.bed > no_aligns_in_house_lowcmp.txt
awk '{sum += $3-$2-1} END {print sum}' no_aligns_sentieon_lowcmp.bed > no_aligns_sentieon_lowcmp.txt
awk '{sum += $3-$2-1} END {print sum}' no_aligns_in_house_lowmap.bed > no_aligns_in_house_lowmap.txt
awk '{sum += $3-$2-1} END {print sum}' no_aligns_sentieon_lowmap.bed > no_aligns_sentieon_lowmap.txt
awk '{sum += $3-$2-1} END {print sum}' no_aligns_in_house_segdup.bed > no_aligns_in_house_segdup.txt
awk '{sum += $3-$2-1} END {print sum}' no_aligns_sentieon_segdup.bed > no_aligns_sentieon_segdup.txt

# get sum of all regions
bedtools merge -i <( bedtools sort -i <(cat no_aligns_in_house_lowcmp.bed \
                                         no_aligns_in_house_lowmap.bed \
                                         no_aligns_in_house_segdup.bed) \
                        -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai) > no_aligns_in_house_all_merged.bed
awk '{sum += $3-$2-1} END {print sum}' no_aligns_in_house_all_merged.bed > no_aligns_in_house_all_merged.txt

bedtools merge -i <( bedtools sort -i <(cat no_aligns_sentieon_lowcmp.bed \
                                         no_aligns_sentieon_lowmap.bed \
                                         no_aligns_sentieon_segdup.bed) \
                        -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai) > no_aligns_sentieon_all_merged.bed
awk '{sum += $3-$2-1} END {print sum}' no_aligns_sentieon_all_merged.bed > no_aligns_sentieon_all_merged.txt



# Create truth difficult regions
# Intersect truth with various difficult regions beds
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions.bed \
    -b Bed_Files/LowComplexity/Human_Full_Genome_TRDB_hg19_150331_all_merged_coordonly.bed \
    > diploid_truth_regions_lowcmp.bed
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions.bed \
    -b Bed_Files/SegmentalDuplications/hg19_self_chain_split_both_coordonly.bed \
    > diploid_truth_regions_segdup.bed
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions.bed \
    -b Bed_Files/mappability/human_g1k_v37_gemmap_l100_m0_e0_nonuniq.sort_coordonly.bed \
    > diploid_truth_regions_lowmap.bed

# Intersection of in house coverage and truth difficult regions
bedtools intersect -a ../Unaligned_Truth_Regions/all_in_house_coverage.bed \
    -b diploid_truth_regions_lowcmp.bed \
    > all_in_house_lowcmp_cov.bed
bedtools intersect -a ../Unaligned_Truth_Regions/all_in_house_coverage.bed \
    -b diploid_truth_regions_segdup.bed \
    > all_in_house_segdup_cov.bed
bedtools intersect -a ../Unaligned_Truth_Regions/all_in_house_coverage.bed \
    -b diploid_truth_regions_lowmap.bed \
    > all_in_house_lowmap_cov.bed

# Intersection of sentieon regions and truth difficult regions
bedtools intersect -a ../Unaligned_Truth_Regions/all_sentieon_coverage.bed \
    -b diploid_truth_regions_lowcmp.bed \
    > all_sentieon_lowcmp_cov.bed
bedtools intersect -a ../Unaligned_Truth_Regions/all_sentieon_coverage.bed \
    -b diploid_truth_regions_segdup.bed \
    > all_sentieon_segdup_cov.bed
bedtools intersect -a ../Unaligned_Truth_Regions/all_sentieon_coverage.bed \
    -b diploid_truth_regions_lowmap.bed \
    > all_sentieon_lowmap_cov.bed

# Intersection of in house uncovered truth and truth difficult regions
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed \
    -b diploid_truth_regions_lowcmp.bed \
    > all_in_house_lowcmp_no_cov.bed
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed \
    -b diploid_truth_regions_segdup.bed \
    > all_in_house_segdup_no_cov.bed
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_in_house.bed \
    -b diploid_truth_regions_lowmap.bed \
    > all_in_house_lowmap_no_cov.bed

# Intersection of sentieon uncovered truth and truth difficult regions
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed \
    -b diploid_truth_regions_lowcmp.bed \
    > all_sentieon_lowcmp_no_cov.bed
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed \
    -b diploid_truth_regions_segdup.bed \
    > all_sentieon_segdup_no_cov.bed
bedtools intersect -a ../Unaligned_Truth_Regions/diploid_truth_regions_no_aligns_sentieon.bed \
    -b diploid_truth_regions_lowmap.bed \
    > all_sentieon_lowmap_no_cov.bed

# merge all regions, then sum the coverage
bedtools merge -i <( bedtools sort -i <(cat all_in_house_lowcmp_cov.bed \
                                         all_in_house_lowmap_cov.bed \
                                         all_in_house_segdup_cov.bed) \
                        -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai) > all_in_house_difficult_regions_merged_cov.bed
awk '{sum += $3-$2-1} END {print sum}' all_in_house_difficult_regions_merged_cov.bed > all_in_house_difficult_regions_merged_cov.txt

bedtools merge -i <( bedtools sort -i <(cat all_sentieon_lowcmp_cov.bed \
                                         all_sentieon_lowmap_cov.bed \
                                         all_sentieon_segdup_cov.bed) \
                        -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai) > all_sentieon_difficult_regions_merged_cov.bed
awk '{sum += $3-$2-1} END {print sum}' all_sentieon_difficult_regions_merged_cov.bed > all_sentieon_difficult_regions_merged_cov.txt

bedtools merge -i <( bedtools sort -i <(cat all_in_house_lowcmp_no_cov.bed \
                                         all_in_house_lowmap_no_cov.bed \
                                         all_in_house_segdup_no_cov.bed) \
                        -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai) > all_in_house_difficult_regions_merged_no_cov.bed
awk '{sum += $3-$2-1} END {print sum}' all_in_house_difficult_regions_merged_no_cov.bed > all_in_house_difficult_regions_merged_no_cov.txt

bedtools merge -i <( bedtools sort -i <(cat all_sentieon_lowcmp_no_cov.bed \
                                         all_sentieon_lowmap_no_cov.bed \
                                         all_sentieon_segdup_no_cov.bed) \
                        -g /research/rv-02/home/eanderson/Resources_And_DBs/hs37d5.fa.fai) > all_sentieon_difficult_regions_merged_no_cov.bed
awk '{sum += $3-$2-1} END {print sum}' all_sentieon_difficult_regions_merged_no_cov.bed > all_sentieon_difficult_regions_merged_no_cov.txt
