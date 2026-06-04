# Alignments

1000 Genomes Project CRAM URL lists used for alignment-based benchmarking.

## Files

- `1000G.hgsvcv3.txt` - 64 CRAM URLs for 1000G samples matching HGSVCv3 assemblies
- `1000G.hprcy1.txt` - 38 CRAM URLs for 1000G samples matching HPRC Year 1 assemblies

Each file contains FTP URLs to CRAM files from the 1000 Genomes Project hosted at EBI.

## Usage

These read sets were used in the leave-zero-out and leave-all-out experiments.
The original CRAMs are aligned to the GRCh38 decoy reference listed in
`../reference/`; benchmark-specific realignment is described in the manuscript
methods and represented by the prediction tables under `../benchmarking/`.

## Source

Original alignments: [1000 Genomes Project at EBI](ftp://ftp.sra.ebi.ac.uk/vol1/run/)
