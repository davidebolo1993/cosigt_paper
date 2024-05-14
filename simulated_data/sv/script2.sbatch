#!/bin/bash
#SBATCH --job-name=wfmash_map
#SBATCH --mail-type=END
#SBATCH --mail-user=davide.bolognini@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=00:25:00
#SBATCH --ntasks=1
#SBATCH --output=%x_%A.log #default to <jobname>_<jobid>.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=30000M

#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

#If you need a conda env 
echo "init conda"
eval "$(conda shell.bash hook)"
conda activate wfmashenv_0130

mkdir -p delpaf
#wfmash ../snp/chr20.fa all_plus_del.fa -X -t 16 > delpaf/all_plus_del.paf
#wfmash ../snp/chr20.fa all_plus_del.fa -X -t 16 -s 10k -p 90 > delpaf/all_plus_del.s10p90.paf
wfmash ../snp/chr20.fa all_plus_del.fa -X -t 16 -s 10k -p 95 > delpaf/all_plus_del.s10p95.paf

echo "Ended: $(date)"
echo "========================
