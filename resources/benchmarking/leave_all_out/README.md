# Leave-All-Out

Prediction results for samples without matching haplotypes in the pangenome graph.

## Structure

- **cmrgs_predictions/** - HPRC Year 1 samples genotyped on HGSVCv3-based graphs
- **svs_predictions/** - HGSVCv3 samples genotyped on HPRC Year 1-based graphs

## Files

- `cosigt_cmrgs_hprcy1_samples_on_hgsvcv3_assemblies.tsv.gz` - HPRC samples on CMRG regions
- `cosigt_svs_hgsvcv3_samples_on_hprcy1_assemblies.tsv.gz` - HGSVC samples on SV regions

## File Format

Tab-separated files with the following columns:

- `sample`, `region`, `gene_name` - Sample and region identifiers
- `hap_1_pred`, `hap_2_pred` - COSIGT predicted haplotypes (best match from alternate pangenome)
- `QV_1_pred`, `QV_2_pred`, `QV_sum_pred` - Quality values for predictions
- `hap_1_best`, `hap_2_best` - Best possible haplotypes from alternate pangenome (ground truth proxy)
- `QV_1_best`, `QV_2_best`, `QV_sum_best` - Quality values for best matches
- `error_rate_*_pred` - Error rates for predictions
- `error_rate_*_best` - Error rates for best matches (theoretical lower bound)
- `error_rate_diff_pred_minus_best` - Difference between predicted and best error rates
- `QV_diff_pred_minus_best` - Difference between predicted and best QV

Haplotype predictions are in the format: `SAMPLE#HAPLOTYPE#CONTIG:START-END`
