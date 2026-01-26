# Benchmarking

Truth sets, benchmarking regions, and COSIGT/Locityper prediction results.

## Structure

- **[cmrgs/](cmrgs)** - Clinically-relevant medically-actionable gene regions
- **[svs/](svs)** - Structural variant regions
- **[flagger/](flagger)** - Regions excluded from benchmarking due to assembly errors
- **[hla_typing/](hla_typing)** - HLA typing comparisons and regions
- **[leave_zero_out/](leave_zero_out)** - Leave-zero-out prediction results (COSIGT and Locityper)
- **[leave_all_out/](leave_all_out)** - Leave-all-out prediction results (COSIGT)

Each subdirectory contains a dedicated `README.md` with detailed information.

## Benchmarking Strategies

### Leave-zero-out

Genotyping samples with matching haplotypes in the pangenome graph at various coverages (1X, 2X, 5X, 30X) and conditions (modern DNA, ancient DNA with and without contamination).

### Leave-all-out

Genotyping samples without matching haplotypes in the pangenome graph using full-coverage sequencing data to assess cross-population performance.
