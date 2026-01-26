# Assemblies

Assembly identifiers and URLs for HGSVC and HPRC pangenome samples used in benchmarking.

## Files

### Assembly URLs

- `hgsvcv3.txt` - URLs to HGSVCv3 phased assembly files (FASTA format)
- `hprcy1.txt` - URLs to HPRC Year 1 phased assembly files (FASTA format)
- `chm13.txt` - URL to CHM13 reference assembly

### Chromosome Assignment

- `hgsvcv3.chroms.txt` - URLs to chromosome assignment files for HGSVCv3 assemblies
- `hprcy1.chroms.txt` - URLs to chromosome assignment files for HPRC assemblies

### Per-Chromosome Assembly Lists

- `hgsvcv3.bychrom/` - Assembly contig identifiers organized by chromosome (chr1-22, X, Y)
- `hprcy1.bychrom/` - Assembly contig identifiers organized by chromosome (chr1-22, X, Y)

Each `chr*.ids.txt` file contains assembly contig names in the [PanSN-spec](https://github.com/pangenome/PanSN-spec) `SAMPLE#HAPLOTYPE#CONTIG` that map to the corresponding reference chromosome.

## Sources

- **HGSVCv3**: [1000 Genomes HGSVC3 collection](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/)
- **HPRC Year 1**: [Human Pangenome Reference Consortium](https://human-pangenomics.s3.amazonaws.com/)
