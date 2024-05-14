#!/bin/bash
module load samtools/1.18
cat ../snp/samples.txt | paste - - - - -  > blocks.txt
seq 1 10 > dim.txt
paste blocks.txt dim.txt > process.txt
rm dim.txt blocks.txt
samtools faidx ../snp/GRCh38_full_analysis_set_plus_decoy_hla.fa chr20 > chr20.fa
bwa-mem2 index chr20.fa
sentence=$(sbatch script1.sbatch)
#wait until all the job ids are consumed
stringarray=($sentence)
jobid=(${stringarray[3]})
VAR=0
while [ "$VAR" -eq 0 ]; do
    sleep 1
    nrow=$(squeue -j $jobid | wc -l)
    if [ $nrow -eq 1 ] ; then #only the header of the jobid table
        VAR=1
    fi
done
cat deletion/samples/*/hap/*fa ../snp/chr20.fa > all_plus_del.fa
sentence=$(sbatch script2.sbatch)
#wait until all the job ids are consumed
stringarray=($sentence)
jobid=(${stringarray[3]})
VAR=0
while [ "$VAR" -eq 0 ]; do
    sleep 1
    nrow=$(squeue -j $jobid | wc -l)
    if [ $nrow -eq 1 ] ; then #only the header of the jobid table
        VAR=1
    fi
done
