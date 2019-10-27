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

      drname="${1/\//}" # strip trailing forward slash
      id="$2"
      t=$3
      mkdir -p paird
      mkdir -p unpaired
      pdr="paired"
      udr="unpaired"
      if [[ $# == 3 ]]; then
         n="$((50/$t))"
         cat $id | sed 's/=/ /g' | awk -v d="$drname" '{print d"/"$1,d"/"$2}' | parallel echo "PE -phred33 {} {} $pdr/{}.fp.gz $pdr/{}.rp.gz $udr/{}.fu.gz $udr/{}.ru.gz ILLUMINACLIP:${HOME}/bioTools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:38 TRAILING:38 SLIDINGWINDOW:4:15 MINLEN:36 -threads $n -o $dr" | xargs -P$n -n17 echo
      else
         cat $id | sed 's/=/ /g' | awk -v d="$drname" '{print d"/"$1,d"/"$2}' | parallel echo "PE -phred33 {} {} $pdr/{}.fp.gz $pdr/{}.rp.gz $udr/{}.fu.gz $udr/{}.ru.gz ILLUMINACLIP:${HOME}/bioTools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:38 TRAILING:38 SLIDINGWINDOW:4:15 MINLEN:36 -threads $n -o $dr" | xargs -P5 -n17 echo
      fi
      echo "Done! All results save in '$dr'"
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
    echo -e """
	Usage: ./trimFastq.sh <dir> <idlist> <threads>

          dir: Path to Fastq files
         list: File containing paired fastq files separated by an equal (=) sign

        e.g. Note the difference in the fastq file names

        ERR1823587_\e[38;5;1m1\e[0m.fastq.gz=ERR1823587_\e[38;5;1m2\e[0m.fastq.gz
        CTRL1_S1_L001_\e[38;5;1mR1\e[0m_001.fastq.gz=CTRL1_S1_L001_\e[38;5;1mR2\e[0m_001.fastq.gz
        ERR1823588_\e[38;5;1m1\e[0m.fastq.gz=ERR1823588_\e[38;5;1m2\e[0m.fastq.gz
    """

fi
