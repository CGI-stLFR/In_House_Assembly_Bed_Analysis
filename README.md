# Bed Analysis for Assemblies

This pipeline is set up to analyse phased or unphased de novo assemblies.
It outputs metrics on unique, multi and total coverage.
It also outputs a heterozygosity summary for missing regions and a summary of stratified uncovered regions.
The stratification summary looks at low complexity, low mappability, and segmental duplication regions.

## Running the pipeline

This pipeline is primarily for simulated libraries.
Create a directory with the name of the sample.
Symlink the fastas as contigs_N.fa, and then add the directory you created to the samples section in the config.
Then run the pipeline as below.

```
# -j supplied threads
snakemake -j 60 2>&1 | tee snakemake.err.txt
```

## Output

Various summaries can be found as detailed below.
Other files contain intermediate output and can potentially be analysed in more depth.

__Genome Coverage Summary__
`Genome_Coverage_Summary/{ref}/{sample}/{hap}/genome_coverage_summary.txt`

__Heterozygosity Summary__
`Heterozygosity_Summary/{ref}/{sample}/heterozygosity_summary.txt`

__Stratification Summary__
`Difficult_Regions_Summary/{ref}/{sample}/summary.txt`

## run_binning.config

- samples
    - samples specifies which samples to evaluate
    - just add the name of the dir that's been added to `Assembly` before running the pipeline
- difficult_regions
    - If you want to analyze other difficult regions, you can modify the dictionary present in the config.
