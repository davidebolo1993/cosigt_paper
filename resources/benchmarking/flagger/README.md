# Flagger

Assembly error masks excluded from benchmarking.

## Files

- `flagger.exclude.hgsvcv3.bed` - HGSVCv3 assembly regions excluded from evaluation
- `flagger.exclude.hprcy2.bed` - HPRC Year 2 assembly regions excluded from evaluation

## Description

Regions flagged by [Flagger](https://github.com/mobinasri/flagger) as assembly
errors, duplications, or collapsed repeats were removed from benchmark
comparisons. This avoids scoring genotyping accuracy in unreliable assembly
contexts.

**Format**: BED file with 3 columns (contig, start, end) using assembly contig names

## Sources

- **HGSVCv3**: [HGSVC3 Flagger results](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20241218_phase3-main-pub_data/uwash/flagger/verkko/final_beds_alt_removed/)
- **HPRC Year 2**: [HPRC intermediate assembly QC](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assembly_qc/flagger)
