
### Data preparation steps

`sbatch` files illustrating input data preparation

- [1](download_verkko.sbatch) - `download_verkko.sbatch`. Download (65) diploid verkko assemblies and adapt the name of the contigs to [PanSN-spec](https://github.com/pangenome/PanSN-spec) format
- [2](download_chm13v2.sh) - `download_chm13v2.sh`. Download chm13v2 reference
- [3](partition_bychrom.sbatch) - `partition_bychrom.sbatch`. Partitionate verkko contigs by chromosome knowing which contig belongs to which chromosome. Also add chm13v2 chromosomes to the per-chromosome verkko contigs.
- [4](download_1000G.sbatch) - `download_1000G.sbatch`. Download short-read samples matching haplotypes in verkko assemblies (short-reads for 64/65 samples). Also get the reference used to decode .cram files and the primary reference for comparison to alternative tools  
- [5](realign_1000G.sbatch) - `realign_1000G.sbatch`. Re-align short-reads to the primary reference
- [6](find_alternative_regions.sh) - `find_alternative_regions.sh`. A small script to include alternative regions in cosigt input bed file given a re-alignment of patch/alternative contigs to the primary reference
