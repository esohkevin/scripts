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
      dname=$2 # strip trailing forward slash
      t=$3
      mkdir -p ${dname}../paired
      mkdir -p ${dname}../unpaired
      pdr="${dname}../paired/"
      udr="${dname}../unpaired/"
      if [[ $# == 3 ]]; then
         n="$((50/$t))"
	 cat $id | parallel echo "PE -phred33 $dname{/}_1.fastq.gz $dname{/}_2.fastq.gz $dname{/}_fp.fastq.gz $dname{/}_fu.fastq.gz $dname{/}_rp.fastq.gz $dname{/}_ru.fastq.gz ILLUMINACLIP:$HOME/bioTools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:38 TRAILING:38 SLIDINGWINDOW:4:15 MINLEN:36 -threads $t" | xargs -P$n -n15 trimmomatic

	 mv ${dname}*_fp.fastq.gz ${dname}*_rp.fastq.gz $pdr
         mv ${dname}*_fu.fastq.gz ${dname}*_ru.fastq.gz $udr

      else
         cat $id | parallel echo "PE -phred33 $dname{/}_1.fastq.gz $dname{/}_2.fastq.gz $dname{/}_fp.fastq.gz $dname{/}_fu.fastq.gz $dname{/}_rp.fastq.gz $dname{/}_ru.fastq.gz ILLUMINACLIP:$HOME/bioTools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:38 TRAILING:38 SLIDINGWINDOW:4:15 MINLEN:36 -threads 1" | xargs -P5 -n15 trimmomatic

	 mv ${dname}*_fp.fastq.gz ${dname}*_rp.fastq.gz $pdr
         mv ${dname}*_fu.fastq.gz ${dname}*_ru.fastq.gz $udr
      fi
      echo "Done! Results saved in '$pdr' and '$udr'"
   fi
else
    echo """
	Usage: ./trimFastq.sh <idlist> <fpath> <threads>

	 idlist: File containing SRA run accessions. (NB: Must be in same path as fastq files)
	  fpath: Path to fastq files. (NB: Must end with a forward '/' slash)
	threads: Optional [default = 1]

	e.g.: File content may be of either format or both
	
	ftp.sra.ebi.ac.uk/vol1/err/ERR182/007/ERR1823587
	ftp.sra.ebi.ac.uk/vol1/err/ERR182/008/ERR1823588
	ftp.sra.ebi.ac.uk/vol1/err/ERR182/009/ERR1823589
	ERR1823590
	ERR1823591
	ERR1823592
    """

fi
