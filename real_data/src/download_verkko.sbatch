#!/bin/bash
#SBATCH --job-name=download.verkko
#SBATCH --mail-type=END
#SBATCH --partition=cpuq
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --cpus-per-task=10
#SBATCH --mem=5G
#SBATCH --array=1-3

set -euo pipefail

#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

#requires bgzip
module load htslib-tools/1.14
module load samtools/1.18

DIR_BASE=$(readlink -f ..)
DATA_BASE=$DIR_BASE/data
RES_BASE=$DIR_BASE/resources
batch=$(head -n $SLURM_ARRAY_TASK_ID $RES_BASE/verkko.batches.txt | tail -1 | cut -f 1)
manifest=$(head -n $SLURM_ARRAY_TASK_ID $RES_BAS/verkko.batches.txt | tail -1 | cut -f 2)
BATCH_BASE=$DATA_BASE/$batch
mkdir -p $BATCH_BASE

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/$batch/$manifest -P $BATCH_BASE

for sample in $(grep ".fasta.gz" $BATCH_BASE/$manifest | grep -v '20230818_verkko_batch1.*HG00732\|HG00732.*20230818_verkko_batch1' | cut -d "/" -f 3 | sort | uniq); do 
    SAMPLE_BASE=$DATA_BASE/$batch/$sample
    mkdir -p $SAMPLE_BASE
    #download to base
    wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/$batch/assemblies/$sample/$sample.vrk-ps-sseq.asm-hap1.fasta.gz -P $SAMPLE_BASE
    wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/$batch/assemblies/$sample/$sample.vrk-ps-sseq.asm-hap2.fasta.gz -P $SAMPLE_BASE
    wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/$batch/assemblies/$sample/$sample.vrk-ps-sseq.asm-unassigned.fasta.gz -P $SAMPLE_BASE
    #to panSN-spec
    zcat $(ls $SAMPLE_BASE/*gz | sort) | sed -e "s/>haplotype1/>$sample#1#haplotype1/g" -e "s/>haplotype2/>$sample#2#haplotype2/g" -e "s/>unassigned/>$sample#U#unassigned/g" | bgzip -@ 10 -l 9 > $SAMPLE_BASE/$sample.fa.gz
    samtools faidx $SAMPLE_BASE/$sample.fa.gz
done

#Closing statement
echo "========================"
echo "Completed: $(date)"
