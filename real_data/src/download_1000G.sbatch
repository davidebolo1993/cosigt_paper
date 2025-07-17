#!/bin/bash
#SBATCH --job-name=download.1000G
#SBATCH --mail-type=END
#SBATCH --partition=cpuq
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G

set -euo pipefail

#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

DIR_BASE=$(readlink -f ..)
DATA_BASE=$DIR_BASE/data
RES_BASE=$DIR_BASE/resources


## READS
for sample in $(grep -h -f <(ls -d $DATA_BASE/*/* | grep -v -E "tsv|txt" | rev | cut -d "/" -f 1 | rev) $RES_BASE/*index | cut -f 1); do
    id=$(echo $sample | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1)
    mkdir -p $DATA_BASE/1000G/$id
    wget -c $sample -P $DATA_BASE/1000G/$id
    wget -c $sample".crai" -P $DATA_BASE/1000G/$id
done

#after this, we still miss NA24385 and NA21487
#NA24385 Project: PRJEB35491
module load samtools/1.18
mkdir -p $DATA_BASE/1000G/NA24385
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR368/ERR3684866/CT21796_NA24385_NIST_UDI22_HiseqX_19122017_Proband_S1.bam -P $DATA_BASE/1000G/NA24385
mv $DATA_BASE/1000G/NA24385/CT21796_NA24385_NIST_UDI22_HiseqX_19122017_Proband_S1.bam $DATA_BASE/1000G/NA24385/NA24385.bam
samtools index $DATA_BASE/1000G/NA24385/NA24385.bam
#this will be realigned to the same reference of the crams

#NA21487?
#This I don't know where to get short-reads from

## REFERENCE(s)
#Reference for decoding crams
#This can be used as hg38 reference for our analyses
mkdir -p $DATA_BASE/reference
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P $DATA_BASE/reference/
bgzip $DATA_BASE/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx $DATA_BASE/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz

#Add locityper reference - we will re-align everything against this for a fair comparison - this is the link they provided in their docs
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz -P $DATA_BASE/reference/
#this is not bgzip-compressed
gzip -d $DATA_BASE/reference/GRCh38.primary_assembly.genome.fa.gz
bgzip $DATA_BASE/reference/GRCh38.primary_assembly.genome.fa
samtools faidx $DATA_BASE/reference/GRCh38.primary_assembly.genome.fa.gz

#since we will re-align to this, index with bwa
module load bwa/0.7.18
bwa index $DATA_BASE/reference/GRCh38.primary_assembly.genome.fa.gz
#the bam will be re-aligned to the decoy reference as well, so we need to index that as well
bwa index $DATA_BASE/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz

#Closing statement
echo "========================"
echo "Completed: $(date)"
