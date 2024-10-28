#!/bin/bash

#see https://davidebolo1993.github.io/cosigtdoc/setup/setup.html
conda activate smk7324app

DIR_BASE=$(readlink -f ..)
RES_BASE=$DIR_BASE/resources
DATA_BASE=$DIR_BASE/data

roi=$RES_BASE/c4.roi.bed
assemblies=$DATA_BASE/verkko_by_chromosome/chr6/chr6.fa.gz
assembliesfa=$DATA_BASE/verkko_by_chromosome/chr6/chr6.fa
zcat $assemblies > $assembliesfa
samples=$DATA_BASE/1000G
reference=$DATA_BASE/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd cosigt/cosigt_smk
#setup cookiecutter profile to run on slurm
#see https://davidebolo1993.github.io/cosigtdoc/setup/setup.html
template="gh:Snakemake-Profiles/slurm"
cookiecutter \
    --output-dir config \
    $template

#prepare sample map
for s in $(ls $samples/*/*.*am); do cram=$(basename $s) && id=$(basename $(dirname $s)) && echo -e "$cram\t$id" >> sample.map.tsv ; done


#download annotations
wget https://raw.githubusercontent.com/davidebolo1993/cosigt/refs/heads/master/docfiles/chr6.c4.annotation.bed

#prepare input
python workflow/scripts/organize.py \
    -a $samples \
    -r $reference \
    --assemblies $assembliesfa \
    --roi $roi \
    --output $DATA_BASE/c4 \
    --samplemap sample.map.tsv \
    --annotation chr6.c4.annotation.bed

#run
sh snakemake.singularity.profile.run.sh

#benchmark
#SINGULARITY_TMPDIR=/tmp snakemake --profile config/slurm --singularity-args "-B /group/soranzo/davide.bolognini/working/dev/cosigt_paper/real_data/data/1000G,/group/soranzo/davide.bolognini/working/dev/cosigt_paper/real_data/data/reference,/group/soranzo/davide.bolognini/working/dev/cosigt_paper/real_data/data/verkko_by_chromosome/chr6,/localscratch,/group/soranzo/davide.bolognini/working/dev/cosigt_paper/real_data/src/cosigt/cosigt_smk,/group/soranzo/davide.bolognini/working/dev/cosigt_paper/real_data/data/c4 -e" benchmark
