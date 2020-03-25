#!/bin/bash

echo """@HD	VN:1.6	SO:coordinate
@SQ	SN:Pf3D7_01_v3	LN:640851
@SQ	SN:Pf3D7_02_v3	LN:947102
@SQ	SN:Pf3D7_03_v3	LN:1067971
@SQ	SN:Pf3D7_04_v3	LN:1200490
@SQ	SN:Pf3D7_05_v3	LN:1343557
@SQ	SN:Pf3D7_06_v3	LN:1418242
@SQ	SN:Pf3D7_07_v3	LN:1445207
@SQ	SN:Pf3D7_08_v3	LN:1472805
@SQ	SN:Pf3D7_09_v3	LN:1541735
@SQ	SN:Pf3D7_10_v3	LN:1687656
@SQ	SN:Pf3D7_11_v3	LN:2038340
@SQ	SN:Pf3D7_12_v3	LN:2271494
@SQ	SN:Pf3D7_13_v3	LN:2925236
@SQ	SN:Pf3D7_14_v3	LN:3291936
@SQ	SN:Pf3D7_API_v3	LN:34250
@SQ	SN:Pf_M76611	LN:5967
@RG	ID:${1}	PL:ILLUMINA	SM:${1}
@PG	ID:bwa	PN:bwa	VN:0.7.12-r1039	CL:bwa mem -t 10 PlasmoDB-39_Pfalciparum3D7_Genome.fasta paired/${1}_fp.fastq.gz paired/${1}_rp.fastq.gz
@PG	ID:MarkDuplicates	VN:2.20.6-SNAPSHOT	CL:MarkDuplicates INPUT=[/BIODATA/StudentsHome/catherineb/genomic_scan/aligned/${1}.bam] OUTPUT=${1}.bam_mkdups.bam METRICS_FILE=${1}.bam_marked_dup_metrics.txt REMOVE_DUPLICATES=true    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag CLEAR_DT=true DUPLEX_UMI=false ADD_PG_TAG_TO_READS=true ASSUME_SORTED=false DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false	PN:MarkDuplicates""" > ${1}.hdr

for i in *bam_mkdups.bam; do echo ${i/.bam_mkdups.bam/.hdr} ${i} ${i/bam_/}; done | parallel --col-sep ' ' echo "reheader -P {1} {2} \\> {3}" | xargs -P50 -n6 echo samtools | parallel

for i in *bam_mkdups.bam; do echo ${i} ${i/bam_/}; done | parallel --col-sep ' ' mv {2} {1}
