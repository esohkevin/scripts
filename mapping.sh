#!/bin/bash

export PATH=${PATH}:/BIODATA/StudentsHome/bin

if [[ $# == 2 ]]; then

    ref="$1"
    paired_path="$2"
    
    fastqBase=`for i in ${paired_path}*forward_paired.fastq.gz ${paired_path}*reverse_paired.fastq.gz; do echo ${i/_*}; done | sort | uniq`
    
    echo -e "\n#~@~# Indexing  Reference sequence #~@~#"
    
    hisat2-build $ref ${ref/.f*/}
    
    echo -e "\n#~@~# Starting alignment #~@~#"
    
       for seqRead in ${fastqBase}; do
    
           hisat2 \
	         -x ${ref/.f*/} \
	         -1 ${seqRead}_forward_paired.fastq.gz \
	         -2 ${seqRead}_reverse_paired.fastq.gz \
	         -S ${seqRead}.sam
    
       done
    
    echo -e "\nAll alignmnets completed: `date` "

else
    echo """
	Usage: ./mapping.sh <ref_seq> <path_to_trimmed_seqs>

    """
fi
