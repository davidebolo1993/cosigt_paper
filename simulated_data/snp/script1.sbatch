#!/bin/bash
#SBATCH --job-name=simulate_snp
#SBATCH --mail-type=END
#SBATCH --mail-user=davide.bolognini@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --output=%x_%A.log #default to <jobname>_<jobid>.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=10000M
#SBATCH --array=1-50


#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

#Load modules
module load bcftools/1.18
module load wgsim/0.0.1
module load bwa-mem2/2.2.1
module load mosdepth/0.3.6
module load bedtools/2.30.0

#If you need a conda env 
echo "init conda"
eval "$(conda shell.bash hook)"
conda activate visorenv_latest

line=$(head -n $SLURM_ARRAY_TASK_ID samples.txt | tail -1)
dim=$(cut -f 2 chr20.fa.fai)
rlen=150
cov=30
mkdir -p samples/$line/hap
mkdir -p samples/$line/fastq
mkdir -p samples/$line/bam
bcftools view --threads 4 -O b -o samples/$line/$line.bcf -s $line -m2 -M2 -c 1 -C 1 ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' samples/$line/$line.bcf | awk -v var=$line 'OFS=FS="\t"''{if($5=="1|0") {print $1, ($2 -1), $2, "SNP", $4, "0" > "samples/"var"/"var".h1.bed"} else if($5 == "0|1"){print $1, ($2 -1), $2, "SNP", $4, "0" > "samples/"var"/"var".h2.bed"}}'
#liftOver samples/$line/$line.h1.bed grch38-chm13v2.chain samples/$line/$line.h1.lifted.bed samples/$line/$line.h1.unmapped.txt
#liftOver samples/$line/$line.h2.bed grch38-chm13v2.chain samples/$line/$line.h2.lifted.bed samples/$line/$line.h2.unmapped.txt
bedtools sort -i samples/$line/$line.h1.bed > samples/$line/$line.h1.sorted.bed
bedtools sort -i samples/$line/$line.h2.bed > samples/$line/$line.h2.sorted.bed
VISOR hack -g chr20.fa -b samples/$line/$line.h1.sorted.bed samples/$line/$line.h2.sorted.bed -o samples/$line/hap
nreads=$(($cov*($dim)/$rlen/2/2)) #per-haplotype number of reads
wgsim samples/$line/hap/h1.fa samples/$line/fastq/r1.h1.fq samples/$line/fastq/r2.h1.fq -1 $rlen -2 $rlen -N $nreads
wgsim samples/$line/hap/h2.fa samples/$line/fastq/r1.h2.fq samples/$line/fastq/r2.h2.fq -1 $rlen -2 $rlen -N $nreads
R1=$(ls samples/$line/fastq/r1* | sort)
R2=$(ls samples/$line/fastq/r2* | sort)
bwa-mem2 mem -t 4 chr20.fa <(cat $R1) <(cat $R2) | samtools sort -@ 4 -o samples/$line/bam/$line.bam --write-index -
mosdepth -n -x --by 500 samples/$line/bam/$line samples/$line/bam/$line.bam
cut -f 4 samples/$line/bam/$line.mosdepth.summary.txt | tail -1 > samples/$line/bam/$line.mean_cov.txt
#rename haplotypes before using them to re-build graph
rm samples/$line/hap/*fai
sed -i "s/>chr20/>$line#1#chr20/" samples/$line/hap/h1.fa && samtools faidx samples/$line/hap/h1.fa
sed -i "s/>chr20/>$line#2#chr20/" samples/$line/hap/h2.fa && samtools faidx samples/$line/hap/h2.fa

#Closing statement
echo "========================"
echo "Completed: $(date)"

