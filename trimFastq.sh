#!/bin/bash


if [[ $# == 1 ]]; then

    #path_to_trim="$1"
    path_to_fq="$1"
    
    fastqBase=$(basename -a $path_to_fq | sort | uniq)
    #trimmomatic="${path_to_trim}trimmomatic-0.39.jar"

    which trimmomatic
    if [[ $? != 0 ]]; then
       echo "Could not locate trimmomatic"
       echo "Please install by one of the following means"
       echo "'sudo apt-get install trimmomatic'"
       echo "'conda install trimmomatic -c bioconda'"
    else

       echo -e """\n\e[38;5;35m#------- Running Trimmomatic -------#\e[0m\n""" 
       
       for seqRead in ${fastqBase}; do
           trimmomatic \
               PE \
               -phred33 \
               ${seqRead} ${seqRead} \
               ${seqRead/fastq.gz/fp.fastq.gz} ${seqRead/fastq.gz/rp.fastq.gz} \
               ${seqRead}_reverse_paired.fastq.gz ${seqRead}_reverse_unpaired.fastq.gz \
               ILLUMINACLIP:${path_to_trim}adapters/TruSeq3-PE.fa:2:30:10 \
               LEADING:3 \
               TRAILING:3 \
               SLIDINGWINDOW:4:15 \
               MINLEN:36 \
               -threads 10
       done
    fi
else
    echo """
	Usage: ./trimFastq.sh <path_to_trommomatic> <path_to_fastq_files>

	NB: All paths should end with a forward slash '/'
	
	To specify the current working directory, enter './'
    """

fi
