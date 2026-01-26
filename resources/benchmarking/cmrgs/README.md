# CMRGs

Clinically-relevant medically-actionable gene regions for benchmarking.

## Files

- `cmrgs.hprcy1.bed` - Original CMRG coordinates (326 genes) as retrieved from the Locityper paper and database
- `cmrgs_refined.hprcy1.bed` - Refined CMRG coordinates after COSIGT refine step (326 genes)

## Description

Both files contain the same 326 genes but with adjusted coordinates. The refined version has coordinates optimized by COSIGT's refinement process to better align with pangenome graph structure.

**Format**: BED file with 4 columns (chromosome, start, end, gene name)