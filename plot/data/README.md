# Plot Data

Input data used by plotting scripts when the data are not already stored in
`../../resources/`.

## Contents

### cyp2d6/

Graph-visualization data for Figure 1 panels C and D:

- `chr22_42031864_42245566.fasta.gz.fai` - Haplotypes in the region
- `chr22_42031864_42245566.clusters.json` - Graph clustering information
- `chr22_42031864_42245566.clusters.medoids.tsv` - Medoid assignments for each cluster
- `chr22_42031864_42245566.node.length.tsv` - Node length information
- `chr22_42031864_42245566.tsv.gz` - Node coverage, with haplotype coverage over graph nodes
- `chr22_42031864_42245566.plot.bed.gz` - Gene coordinates on each haplotype

These files are generated during COSIGT graph visualization runs.

### hla/

HLA typing summary tables for the HLA comparison panels:

- `hla_comparison.tsv` - HLA comparison summary across seven loci
- `hla_comparison_filt.tsv` - Filtered summary excluding T1K alleles with quality <= 0

**Columns:**

- `gene`, `region` - HLA gene and genomic region
- `samples_with_data` - Samples with successful predictions from COSIGT, T1K, and IMP:02
- `samples_with_data_poss` - Subset where the HLA type has at least one representative in the graph
- `sample_accuracy_*`, `haplotype_accuracy_*` - Accuracy metrics computed on `samples_with_data`
- `sample_accuracy_poss_*`, `haplotype_accuracy_poss_*` - Accuracy metrics computed on `samples_with_data_poss`
- Suffix `_cosigt` or `_t1k` indicates the method

### time_mem_benchmark/

Runtime and memory summaries used by the pre-generated runtime/resource panels:

- `benchmark_rule_summary.tsv` - Per-rule benchmark statistics
- `critical_path_analysis.tsv` - Workflow critical path analysis referenced in the manuscript runtime section

**benchmark_rule_summary.tsv columns:**

- `rule`, `rule_clean` - Snakemake rule name and cleaned version
- `n_jobs` - Number of job executions
- `time_*_min` - Runtime statistics in minutes
- `time_min_job`, `time_max_job` - Job IDs with minimum and maximum runtime
- `mem_*_mb` - Memory usage statistics in MB
- `mem_min_job`, `mem_max_job` - Job IDs with minimum and maximum memory usage

### region_distances/

Region-diversity data for the QV-vs-diversity supplementary panel:

- `region_distances_CMRGs.tsv` - Average normalized distance between clusters for CMRGs across tested conditions
