# Leave-Zero-Out

Prediction results for samples with matching haplotypes in the pangenome graph.

## Structure

- **cmrgs_predictions/** - Predictions on CMRG regions (HPRC Year 1 samples)
  - `modern/` - Modern DNA at various coverages (1X, 2X, 5X, 30X)
  - `aDNA/` - Ancient DNA simulations at 1X and 2X with 0% and 10% contamination
- **svs_predictions/** - Predictions on SV regions (HGSVC samples at 30X)

## File Naming

Files follow the pattern: `cosigt_locityper_{region}_{pangenome}_{coverage}[_{condition}].tsv.gz`

- `region`: cmrgs or svs
- `pangenome`: hprcy1 or hgsvcv3 (or subset for aDNA)
- `coverage`: 1X, 2X, 5X, or 30X
- `condition`: (aDNA only) cont0pct or cont10pct

## File Format

Tab-separated files with the following columns:

- `sample`, `region`, `gene_name` - Sample and region identifiers
- `hap_1_pred_locityper`, `hap_2_pred_locityper` - Locityper haplotype predictions
- `hap_1_pred_cosigt`, `hap_2_pred_cosigt` - COSIGT haplotype predictions
- `QV_1_*`, `QV_2_*`, `QV_sum_*` - Quality values (Phred-scaled) for each tool
- `error_rate_1_*`, `error_rate_2_*`, `avg_error_rate_*` - Per-haplotype and average error rates
- `error_rate_diff_cosigt_minus_locityper` - Error rate difference (COSIGT - Locityper)
- `QV_diff_cosigt_minus_locityper` - QV difference (COSIGT - Locityper)

Haplotype predictions are in the format: `SAMPLE#HAPLOTYPE#CONTIG:START-END`

## Usage

These files contain genotyping predictions and accuracy metrics used to generate figures comparing COSIGT and Locityper performance across different coverage levels and conditions.
