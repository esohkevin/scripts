#!/bin/bash

#### gatk CreateSequenceDictionary -R ref.fasta #### How to create fasta index with gatk when it is needed
re="/BIODATA/StudentsHome/catherineb/genomic_scan/PlasmoDB-39_Pfalciparum3D7_Genome.fasta"
bam="/BIODATA/StudentsHome/catherineb/genomic_scan/aligned/"
ks="/BIODATA/StudentsHome/catherineb/genomic_scan/3d7_hb3.combined.final.vcf"

conda activate cathyinformatics

gatk CreateSequenceDictionary -R $re
#for i in cam sa; do cat ${i}.bamlist.fof | parallel echo --threads 15 -o ${i}fqout {} | xargs -P5 -n5 fastqc; done

sed 's/.bam//g' bqsr.input.list | parallel echo MarkDuplicates I=${bam}{}.bam O={}_mkdups.bam M={}_marked_dup_metrics.txt REMOVE_DUPLICATES=true | xargs -P5 -n5 picard
sed 's/.bam//g' bqsr.input.list | parallel echo bgzip {}_mkdups.bam
sed 's/.bam//g' bqsr.input.list | parallel echo -T RealignerTargetCreator -R $re --known $ks -I {}_mkdups.bam -o {}.intervals | xargs -P10 -n10 gatk3
sed 's/.bam//g' bqsr.input.list | parallel echo -T IndelRealigner -R $re -known $ks -I {}_mkdups.bam -targetIntervals {}.intervals -o {}.realigned.bam | xargs -P10 -n12 gatk3
sed 's/.bam//g' bqsr.input.list | parallel echo bgzip {}.realigned.bam
sed 's/.bam//g' bqsr.input.list | parallel echo BaseRecalibrator -I {}.realigned.bam -R $re --known-sites $ks -O {}_recal_data.table | xargs -P5 -n9 gatk
sed 's/.bam//g' bqsr.input.list | parallel echo ApplyBQSR -R $re -I {}.realigned.bam --bqsr-recal-file {}_recal_data.table -O {}.bam | xargs -P5 -n9 gatk

#cat bqsr.input.list | parallel echo  HaplotypeCaller -I {} -R $re --dbsnp $dbsnp2 -RF MappingQualityAvailableReadFilter -RF MappingQualityNotZeroReadFilter -RF MappedReadFilter -RF NotDuplicateReadFilter -RF NotSecondaryAlignmentReadFilter -RF NonZeroFragmentLengthReadFilter -O {}_gatk.sing.call.vcf.gz | xargs -P10 -n21 gatk
#gatk MergeVcfs -I merge.vcfs.list -O test.vcf.gz

#gatk HaplotypeCaller -I cam.bamlist.list -R $re --dbsnp $dbsnp2 -RF MappingQualityAvailableReadFilter -RF MappingQualityNotZeroReadFilter -RF MappedReadFilter -RF NotDuplicateReadFilter -RF NotSecondaryAlignmentReadFilter -RF NonZeroFragmentLengthReadFilter --minimum-mapping-quality 2 -O cam_gatk.vcf.gz
#####new
