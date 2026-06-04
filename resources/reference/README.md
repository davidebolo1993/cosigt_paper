# Reference

GRCh38 reference genome URL manifests used for alignment, realignment, and
benchmarking.

## Files

- `grch38.decoy.txt` - One URL for the GRCh38 full analysis set with decoys and HLA contigs
- `grch38.primary.txt` - One URL for the GRCh38 primary assembly from GENCODE release 46

## Usage

- `grch38.decoy.txt` records the original reference used by the 1000 Genomes CRAM alignments listed in `../alignments/`.
- `grch38.primary.txt` records the primary assembly used when samples were realigned for the benchmarking experiments.

## Sources

- **Decoy reference**: [1000 Genomes GRCh38 full analysis set](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)
- **Primary assembly**: [GENCODE release 46](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/)
