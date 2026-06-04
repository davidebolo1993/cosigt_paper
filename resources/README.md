# Resources

Resource manifests, benchmarking regions, assembly masks, and prediction results
used in the COSIGT manuscript.

This directory is organized around the inputs and outputs needed to understand
the benchmark tables used by the plotting scripts. Raw read and assembly files
are not stored here; URL manifests point to the external sources.

## Structure

- **[alignments/](alignments)** - 1000 Genomes CRAM URL lists for samples matching HGSVCv3 and HPRC Year 1 assemblies
- **[assemblies/](assemblies)** - HGSVCv3, HPRC Year 1, CHM13, and per-chromosome assembly manifests
- **[benchmarking/](benchmarking)** - Benchmark regions, exclusion masks, HLA resources, and COSIGT/Locityper prediction tables
- **[reference/](reference)** - GRCh38 reference URL manifests

Each subdirectory contains a dedicated `README.md` with more detailed information.
