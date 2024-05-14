#!/usr/bin/bash
module load bcftools/1.18
module load bwa-mem2/2.2.1
module load samtools/1.18
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr20 > chr20.fa
samtools faidx chr20.fa
#keep 50 random samples
bcftools query -l ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz | shuf --random-source chr20.fa | head -50 > samples.txt
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
sed -i "s/^>chr20/>grch38#chr20/" chr20.fa && samtools faidx chr20.fa
cat chr20.fa samples/*/hap/*.fa > all.fa && samtools faidx all.fa
sentence=$(sbatch script3.sbatch)
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
conda activate renv_latest
Rscript snpfreq.r
