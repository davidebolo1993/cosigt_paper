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

#Closing statement
echo "========================"
echo "Completed: $(date)"