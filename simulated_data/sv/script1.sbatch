#!/bin/bash
#SBATCH --job-name=simulate_snp
#SBATCH --mail-type=END
#SBATCH --mail-user=davide.bolognini@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --output=%x_%A.log #default to <jobname>_<jobid>.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=10000M
#SBATCH --array=1-10

#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

#Load modules
module load wgsim/0.0.1
module load bwa-mem2/2.2.1
module load mosdepth/0.3.6
module load samtools/1.18

#If you need a conda env 
echo "init conda"
eval "$(conda shell.bash hook)"
conda activate visorenv_latest

line=$(head -n $SLURM_ARRAY_TASK_ID process.txt | tail -1)
rlen=150
cov=30
dim=$(echo $line | cut -d " " -f 6 | sed 's/$/000/')
end=$((400000 + $dim))
sample=$(echo $line | cut -d " " -f 1-5)
#hap 1 for this samples will have a $dim-dimension alteration, hap2 will be the same

for s in $sample; do

    mkdir -p deletion/samples/$s/hap
    mkdir -p deletion/samples/$s/fastq
    mkdir -p deletion/samples/$s/bam
    echo -e "$s#1#chr20\t400000\t$end\tdeletion\tNone\t0\n" >  deletion/samples/$s/h1.bed
    VISOR hack -g ../snp/samples/$s/hap/h1.fa -b deletion/samples/$s/h1.bed -o deletion/samples/$s/hap
    cp ../snp/samples/$s/hap/h2.fa* deletion/samples/$s/hap/.
    nb1=$(cut -f 2 deletion/samples/$s/hap/h1.fa.fai)
    nb2=$(cut -f 2 deletion/samples/$s/hap/h2.fa.fai)
    nreads1=$(($cov*($nb1)/$rlen/2/2))
    nreads2=$(($cov*($nb2)/$rlen/2/2))
    wgsim deletion/samples/$s/hap/h1.fa deletion/samples/$s/fastq/r1.h1.fq deletion/samples/$s/fastq/r2.h1.fq -1 $rlen -2 $rlen -N $nreads1
    wgsim deletion/samples/$s/hap/h2.fa deletion/samples/$s/fastq/r1.h2.fq deletion/samples/$s/fastq/r2.h2.fq -1 $rlen -2 $rlen -N $nreads2
    R1=$(ls deletion/samples/$s/fastq/r1* | sort)
    R2=$(ls deletion/samples/$s/fastq/r2* | sort)
    bwa-mem2 mem -t 4 chr20.fa <(cat $R1) <(cat $R2) | samtools sort -@ 4 -o deletion/samples/$s/bam/$s.bam --write-index -
    mosdepth -n -x --by 500 deletion/samples/$s/bam/$s deletion/samples/$s/bam/$line.bam
    cut -f 4 deletion/samples/$s/bam/$s.mosdepth.summary.txt | tail -1 > deletion/samples/$s/bam/$s.mean_cov.txt               

done

#Closing statement
echo "========================"
echo "Completed: $(date)"
