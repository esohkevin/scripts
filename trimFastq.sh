#!/bin/bash

if [[ $# == [23] ]]; then

   which trimmomatic
   if [[ $? != 0 ]]; then
      echo "Could not locate trimmomatic"
      echo "Please install by one of the following means"
      echo "'sudo apt-get install trimmomatic'"
      echo "'conda install trimmomatic -c bioconda'"
      echo "Add bioconda channel if it doesn't already exist: 'conda config --add chennels bioconda'"
   else
      id=$1
      dname="$(dirname $1)" # strip trailing forward slash
      s="$2"
      t=$3
      mkdir -p paird
      mkdir -p unpaired
      pdr="paired"
      udr="unpaired"
      if [[ $# == 3 ]]; then
         n="$((50/$t))"
         cat $id | parallel --colsep "$s" echo "PE -phred33 $dname/{1} $dname/{2} $pdr/{1.}.fp.gz $pdr/{2.}.rp.gz $udr/{1.}.fu.gz $udr/{2.}.ru.gz ILLUMINACLIP:${HOME}/bioTools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:38 TRAILING:38 SLIDINGWINDOW:4:15 MINLEN:36 -threads $n" | xargs -P$n -n15 trimmomatic
      else
         cat $id | parallel --colsep "$s" echo "PE -phred33 $dname/{1} $dname/{2} $pdr/{1.}.fp.gz $pdr/{2.}.rp.gz $udr/{1.}.fu.gz $udr/{2.}.ru.gz ILLUMINACLIP:${HOME}/bioTools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:38 TRAILING:38 SLIDINGWINDOW:4:15 MINLEN:36 -threads 1" | xargs -P5 -n15 trimmomatic

      fi
      echo "Done! Results saved in '$pdr' and '$udr'"
   fi
# if [[ $# == 1 ]]; then
# 
#     #path_to_trim="$1"
#     path_to_fq="$1"
#     
#     fastqBase=$(basename -a $path_to_fq | sort | uniq)
#     #trimmomatic="${path_to_trim}trimmomatic-0.39.jar"
# 
#     which trimmomatic
#     if [[ $? != 0 ]]; then
#        echo "Could not locate trimmomatic"
#        echo "Please install by one of the following means"
#        echo "'sudo apt-get install trimmomatic'"
#        echo "'conda install trimmomatic -c bioconda'"
#     else
# 
#        echo -e """\n\e[38;5;35m#------- Running Trimmomatic -------#\e[0m\n""" 
#        
#        for seqRead in ${fastqBase}; do
#            trimmomatic PE -phred33 ${seqRead} ${seqRead} ${seqRead/fastq.gz/fp.fastq.gz} ${seqRead/fastq.gz/rp.fastq.gz} ${seqRead}_reverse_paired.fastq.gz ${seqRead}_reverse_unpaired.fastq.gz ILLUMINACLIP:${path_to_trim}adapters/TruSeq3-PE.fa:2:30:10 LEADING:38 TRAILING:38 SLIDINGWINDOW:4:15 MINLEN:36 -threads 10
#        done
#     fi
else
    echo """
	Usage: ./trimFastq.sh <idlist> <sep> <threads>

        idlist: File containing paired fastq files. Should be in the same path as the Fastq files
           sep: idlist separator

        e.g. Note the difference in the fastq file names

	ERR1823587_1.fastq.gz	ERR1823587_2.fastq.gz : tab ('\t') separated 
	CTRL1_S1_L001_R1_001.fastq.gz;CTRL1_S1_L001_R2_001.fastq.gz : semi-column (';') separated
	ERR1823588_1.fastq.gz ERR1823588_2.fastq.gz : space (' ') separated

	Example:
	trimFastq.sh list.txt '\t' 5
    """

fi
