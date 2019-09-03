#!/bin/bash


if [[ $# == 2 ]]; then

    path_to_trim="$1"
    path_to_fq="$2"
    
    fastqBase=`for i in ${path_to_fq}*1.fastq.gz ${path_to_fq}*2.fastq.gz; do echo ${i/_*}; done | sort | uniq`
    trimmomatic="${path_to_trim}trimmomatic-0.39.jar"
    
    echo -e """\n\e[38;5;35m#------- Running Trimmomatic -------#\e[0m\n""" 
    
    for seqRead in ${fastqBase}; do
        java -jar $trimmomatic \
            PE \
            -phred33 \
            ${seqRead}_1.fastq.gz ${seqRead}_2.fastq.gz \
            ${seqRead}_forward_paired.fastq.gz ${seqRead}_forward_unpaired.fastq.gz \
            ${seqRead}_reverse_paired.fastq.gz ${seqRead}_reverse_unpaired.fastq.gz \
            ILLUMINACLIP:${path_to_trim}adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            -threads 10
    done

else
    echo """
	Usage: ./trimFastq.sh <path_to_trommomatic> <path_to_fastq_files>

	NB: All paths should end with a forward slash '/'
	
	To specify the current working directory, enter './'
    """

fi
