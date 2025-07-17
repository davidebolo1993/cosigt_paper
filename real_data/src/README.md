
### Data preparation steps

`sbatch` files illustrating input data preparation

- [1](download_verkko.sbatch) - `download_verkko.sbatch`. Download (65) diploid verkko assemblies and adapt the name of the contigs to [PanSN-spec](https://github.com/pangenome/PanSN-spec) format
- [2](partition_bychrom.sbatch) - `partition_bychrom.sbatch`. Partitionate verkko contigs by chromosome knowing which contig belongs to which chromosome
- [3](download_1000G.sbatch) - `download_1000G.sbatch`. Download short-read samples matching haplotypes in verkko assemblies (short-reads for 64/65 samples). Also get the reference used to decode .cram files and the primary reference for comparison to alternative tools  
- [4](realign_1000G.sbatch) - `realign_1000G.sbatch`. Re-align short-reads to the primary reference
