# Plot Data

Complementary input data for generating manuscript figures.

## Contents

### cyp2d6/

Data files for Figure 1 panels C and D:

- `chr22_42031864_42245566.fasta.gz.fai` - Haplotypes in the region
- `chr22_42031864_42245566.clusters.json` - Graph clustering information
- `chr22_42031864_42245566.clusters.medoids.tsv` - Medoid assignments for each cluster
- `chr22_42031864_42245566.node.length.tsv` - Node length information
- `chr22_42031864_42245566.tsv.gz` - Node coverage (haplotype coverage over nodes)
- `chr22_42031864_42245566.plot.bed.gz` - Gene coordinates on each haplotype

These files are routinely generated within every COSIGT run.

### hla/

HLA typing comparison data:

- `hla_comparison.tsv` - Full HLA comparison results (Extended Figure 1)
- `hla_comparison_filt.tsv` - Filtered comparison excluding T1K alleles with quality <= 0 (Supplementary Data Figure 9)

**Columns:**
- `gene`, `region` - HLA gene and genomic region
- `samples_with_data` - Samples with successful predictions from both COSIGT and T1K
- `samples_with_data_poss` - Subset of above where HLA types have at least one representative in the graph (used for Extended Data Figure 1 statistics)
- `sample_accuracy_*`, `haplotype_accuracy_*` - Accuracy metrics computed on `samples_with_data`
- `sample_accuracy_poss_*`, `haplotype_accuracy_poss_*` - Accuracy metrics computed on `samples_with_data_poss`
- Suffix `_cosigt` or `_t1k` indicates the method

### time_mem_benchmark/

Runtime and memory usage summaries generated alongside Supplementary Data Figure 8:

- `benchmark_rule_summary.tsv` - Per-rule benchmark statistics
- `critical_path_analysis.tsv` - Workflow critical path analysis (referenced in manuscript section: `Runtime and memory usage`)

**benchmark_rule_summary.tsv columns:**
- `rule`, `rule_clean` - Snakemake rule name and cleaned version
- `n_jobs` - Number of job executions
- `time_*_min` - Runtime statistics (min, mean, median, max) in minutes
- `time_min_job`, `time_max_job` - Job IDs with minimum and maximum runtime
- `mem_*_mb` - Memory usage statistics (min, mean, median, max) in MB
- `mem_min_job`, `mem_max_job` - Job IDs with minimum and maximum memory usage
