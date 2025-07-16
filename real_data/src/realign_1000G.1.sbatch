#!/bin/bash
#SBATCH --job-name=realign.1000G
#SBATCH --mail-type=END
#SBATCH --partition=cpuq
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --array=1-1

set -euo pipefail

#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

module load samtools/1.18
module load bwa/0.7.18

DIR_BASE=$(readlink -f ..)
DATA_BASE=$DIR_BASE/data
REF_DECOY=$DATA_BASE/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
REF_PRIMARY=$DATA_BASE/reference/GRCh38.primary_assembly.genome.fa.gz
aln=$(ls $DATA_BASE/1000G/*/*am | sort | head -n $SLURM_ARRAY_TASK_ID | tail -1)
sample=$(echo $aln | rev | cut -d "/" -f 2 | rev)
outdir=$DATA_BASE/1000G_realigned/$sample
mkdir -p $outdir
#sort cram by name - necessary if we want to get pairs
samtools sort -@ 15 -n --reference $REF_DECOY -o $outdir/$sample.bam -T $outdir/$sample $aln
#fastq pairs
samtools fastq -@ 15 -1 $outdir/$sample.r1.fq.gz -2 $outdir/$sample.r2.fq.gz -0 $outdir/$sample.other.fq.gz -s $outdir/$sample.se.fq.gz -F 0x0 $outdir/$sample.bam
rm $outdir/$sample.bam
#re-align to primary, sort, compress to cram
#index has been generated before
#bwa mem \
	#-t 15 \
	#-Y -K 100000000 \
	#$REF_PRIMARY \
	#$outdir/$sample.r1.fq.gz \
	#$outdir/$sample.r2.fq.gz | samtools sort -@ 10 -m 5G -T $outdir/$sample | samtools view -@ 10 -T $REF_PRIMARY -C -o $outdir/$sample.cram --write-index

#rm $outdir/$sample.*gz
#Closing statement
echo "========================"
echo "Completed: $(date)"
