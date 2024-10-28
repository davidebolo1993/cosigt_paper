#!/bin/bash
#SBATCH --job-name=download.1000G
#SBATCH --mail-type=END
#SBATCH --partition=cpuq
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --output=%x_%A_%a.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

DIR_BASE=$(readlink -f ..)
DATA_BASE=$DIR_BASE/data
RES_BASE=$DIR_BASE/resources

for sample in $(grep -h -f <(ls -d $DATA_BASE/*/* | grep -v -E "tsv|txt" | rev | cut -d "/" -f 1 | rev) $RES_BASE/*index | cut -f 1); do
    id=$(echo $sample | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1)
    echo $id
    mkdir -p $DATA_BASE/1000G/$id
    wget -c $sample -P $DATA_BASE/1000G/$id
    wget -c $sample".crai" -P $DATA_BASE/1000G/$id
done

#after that, we still miss NA24385 and NA21487

#NA24385 Project: PRJEB35491
module load samtools/1.18
mkdir -p $DATA_BASE/1000G/NA24385
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR368/ERR3684866/CT21796_NA24385_NIST_UDI22_HiseqX_19122017_Proband_S1.bam -P $DATA_BASE/1000G/NA24385
mv $DATA_BASE/1000G/NA24385/CT21796_NA24385_NIST_UDI22_HiseqX_19122017_Proband_S1.bam $DATA_BASE/1000G/NA24385/NA24385.bam
samtools index $DATA_BASE/1000G/NA24385/NA24385.bam

#NA21487?

#also download reference for decoding crams
mkdir -p $DATA_BASE/reference
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P $DATA_BASE/reference/


#Closing statement
echo "========================"
echo "Completed: $(date)"
