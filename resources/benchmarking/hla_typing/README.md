# HLA Typing

HLA typing comparison data and genomic regions.

## Files

- `hla_comparison_all_samples.tsv.gz` - Full HLA typing comparison results (excluded from repository)
- `hla_walt.bed` - HLA gene regions with alternative haplotypes
- `hprcy2.chr6.txt` - HPRC Year 2 chromosome 6 contig identifiers

## Description

### hla_walt.bed

BED file containing HLA gene regions (HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1) with:
- Reference coordinates (chr6) for primary assembly
- Gene name
- Alternative contig coordinates from GRCh38 decoy

**Format**: BED file with 5 columns (chromosome, start, end, gene, allele_list)

### hprcy2.chr6.txt

List of HPRC Year 2 assembly contig identifiers mapping to chromosome 6, used for constructing the HLA-specific pangenome graph.

**Format**: One contig identifier per line in the  [PanSN-spec](https://github.com/pangenome/PanSN-spec) `SAMPLE#HAPLOTYPE#CONTIG`

## Note

`hla_comparison_all_samples.tsv.gz` contains sensitive sample-level data and is excluded from the repository. Summary statistics are available in the plotting data directory.
