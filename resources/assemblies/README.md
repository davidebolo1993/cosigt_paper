# Assemblies

Assembly URL manifests and chromosome-specific contig lists for the pangenomes
used in benchmarking.

## Files

### Assembly URLs

- `hgsvcv3.txt` - 195 HGSVCv3 phased assembly FASTA URLs
- `hprcy1.txt` - 76 HPRC Year 1 phased assembly FASTA URLs
- `chm13.txt` - CHM13 reference assembly URL

### Chromosome Assignment

- `hgsvcv3.chroms.txt` - URLs to chromosome-assignment reports for HGSVCv3 assemblies
- `hprcy1.chroms.txt` - URLs to chromosome-assignment reports for HPRC Year 1 assemblies

### Per-Chromosome Assembly Lists

- `hgsvcv3.bychrom/` - Assembly contig identifiers organized by chromosome (chr1-22, X, Y)
- `hprcy1.bychrom/` - Assembly contig identifiers organized by chromosome (chr1-22)

Each `chr*.ids.txt` file contains assembly contig names in the [PanSN-spec](https://github.com/pangenome/PanSN-spec) `SAMPLE#HAPLOTYPE#CONTIG` that map to the corresponding reference chromosome.

## Sources

- **HGSVCv3**: [1000 Genomes HGSVC3 collection](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/)
- **HPRC Year 1**: [Human Pangenome Reference Consortium](https://human-pangenomics.s3.amazonaws.com/)
