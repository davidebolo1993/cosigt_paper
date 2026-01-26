# Flagger

Assembly error regions excluded from benchmarking.

## Files

- `flagger.exclude.hgsvcv3.bed` - Error regions in HGSVCv3 assemblies
- `flagger.exclude.hprcy2.bed` - Error regions in HPRC Year 2 assemblies

## Description

Regions flagged by [Flagger](https://github.com/mobinasri/flagger) as containing assembly errors (Err), duplications (Dup), or collapsed repeats (Col). These regions were excluded from all benchmarking analyses to avoid evaluating genotyping accuracy in unreliable assembly contexts.

**Format**: BED file with 3 columns (contig, start, end) using assembly contig names

## Sources

- **HGSVCv3**: [HGSVC3 Flagger results](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20241218_phase3-main-pub_data/uwash/flagger/verkko/final_beds_alt_removed/)
- **HPRC Year 2**: [HPRC intermediate assembly QC](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assembly_qc/flagger)

