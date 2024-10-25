#!/bin/bash
#SBATCH --job-name=partition.verkko
#SBATCH --mail-type=END
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --output=%x_%A_%a.log
#SBATCH --cpus-per-task=5
#SBATCH --mem=1G
#SBATCH --array=1-24

set -euo pipefail

#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

module load samtools/1.18
module load htslib-tools/1.14 

DIR_BASE=$(readlink -f ..)
DATA_BASE=$DIR_BASE/data
RES_BASE=$DIR_BASE/resources
chrom=$(seq 1 22 | awk '{print "chr"$1} END {print "chrX\nchrY"}' | head -n $SLURM_ARRAY_TASK_ID | tail -1)
mkdir -p $DATA_BASE/verkko_by_chromosome/$chrom
echo $chrom
for hap in $(zgrep $chrom $RES_BASE/"HGSVC3.partitions.tsv.gz" | cut -f 1); do
    sampleid=$(echo $hap | cut -d "#" -f 1)
    sampledir=$(find $DATA_BASE/*verkko* -name $sampleid)
    echo -e $hap"\t"$sampleid"\t"$sampledir
    samtools faidx $sampledir/$sampleid".fa.gz" $hap >> $DATA_BASE/verkko_by_chromosome/$chrom/$chrom".fa"
done

bgzip -@ 5 -l 9 $DATA_BASE/verkko_by_chromosome/$chrom/$chrom".fa"
samtools faidx $DATA_BASE/verkko_by_chromosome/$chrom/$chrom".fa.gz"