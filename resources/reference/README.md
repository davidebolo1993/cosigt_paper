# Reference

GRCh38 reference genome URLs used for alignment and benchmarking.

## Files

- `grch38.decoy.txt` - GRCh38 full analysis set with decoy sequences and HLA contigs
- `grch38.primary.txt` - GRCh38 primary assembly (GENCODE release 46)

## Usage

- **grch38.decoy.txt**: Original reference genome for 1000 Genomes Project alignments (see `../alignments/`)
- **grch38.primary.txt**: Reference genome used for realigning samples in leave-zero-out and leave-all-out benchmarking experiments

## Sources

- **Decoy reference**: [1000 Genomes GRCh38 full analysis set](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)
- **Primary assembly**: [GENCODE release 46](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/)
