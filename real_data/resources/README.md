### Resources description

A brief description of small files included in this sub-folder.
These files are mainly used to download/organize other resources in the [data](../data) folder.

- [1](1000G_2504_high_coverage.sequence.index) - `1000G_2504_high_coverage.sequence.index`. Download [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index)
- [2](1000G_698_related_high_coverage.sequence.index) - `1000G_698_related_high_coverage.sequence.index`. Download [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index)
- [3](verkko.batches.txt) - `verkko.batches.txt`. List of verkko batches.
- Batch1 manifest is [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20230818_verkko_batch1/20230818_verkko_batch1.MANIFEST.txt).
- Batch2 manifest is [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20230927_verkko_batch2/20230927_verkko_batch2.MANIFEST.tsv).
- Batch3 manifest is [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20240201_verkko_batch3/20240201_verkko_batch3.manifest.tsv).
- [4](HGSVC3.partitions.tsv.gz) - `HGSVC3.partitions.tsv.gz`. HGSVCv3 contigs (header follows PanSN-spec) with the matching reference chromosome - determined by re-aligning each contig to [CHM13v2.0](https://github.com/marbl/CHM13).
- [5](GRCh38.alts_to_primary.paf.gz*) - `GRCh38.alts_to_primary.paf.gz*`. Bgzip-compressed alignment of GRCh38 alternative sequences to the primary chromosomes (and their index) - generated with `minimap2 -x asm20 --eqx -c`. Can be used within `impg` to identify intervals in alternative sequences that may have captured reads in the primary chromosomes. 
