#!/bin/bash

module load samtools/1.18

DIR_BASE=$(readlink -f ..)
DATA_BASE=$DIR_BASE/data
mkdir -p $DATA_BASE/chm13v2
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz -P $DATA_BASE/chm13v2
samtools faidx $DATA_BASE/chm13v2/chm13v2.0_maskedY_rCRS.fa.gz
