#!/bin/bash
#PBS -l select=1:ncpus=24
#PBS -l walltime=96:00:00
#PBS -q smp
#PBS -P CBBI1243
#PBS -o /mnt/lustre/groups/CBBI1243/KEVIN/camhi/stdout.txt
#PBS -e /mnt/lustre/groups/CBBI1243/KEVIN/camhi/stderr.txt
#PBS -N JointCalling
#PBS -M kevin.esoh@students.jkuat.ac.ke
#PBS -m b

indir="/mnt/lustre/groups/1000genomes/annotation/annovar/"
anntable="${indir}table_annovar.pl"
annvar="${indir}annotate_variation.pl"
convann="${indir}convert2annovar.pl"
bgzip="/mnt/lustre/users/kesoh/bioTools/bgzip"
exmdir="${indir}example/"
hdb="${indir}humandb/"
invcf="ex2.vcf"
myvcf="joint_call_new.vcf"
snplist="rsids.txt"
ref="/mnt/lustre/groups/1000genomes/annotation/REF_VC/human_g1k_v37.fasta"
re="/mnt/lustre/groups/CBBI1243/Data/hg19/ucsc.hg19.fasta"

module add chpc/BIOMODULES
module load bcftools/1.6.33 samtools/1.9
module add plink/1.90
module add gnu/parallel-20160422

cd /mnt/lustre/groups/CBBI1243/KEVIN/camhi/

# for i in $(cat bamlist.fof); do 
#    b=$(basename $i); samtools view $i -h -T $ref | sed 's/chr//g' > ${b/.bam/U.sam}; 
#    samtools view -h ${b/.bam/U.sam} -O BAM -o ${b/.bam/U.bam}
#    rm ${b/.bam/U.sam}
# done
#bcftools mpileup --thread 15 -f $re -b bamlist.fof | bcftools call -mv -Oz -o cam.vcf.gz
#for i in cam; do echo $i; done | parallel echo "mpileup --thread 15 -f $re -b {}.bamlist.fof \\| bcftools call -mv -Oz -o {}.vcf.gz" | xargs -P2 -n14 bcftools
bcftools mpileup --min-MQ 2 --thread 30 -f $re -Oz -o cam.vcf.gz -b cam.bamlist.fof
bcftools call -mv --threads 30 -Oz -o cam.joint_call.vcf.gz cam.vcf.gz
