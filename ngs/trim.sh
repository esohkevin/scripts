#!/bin/bash

#-----------------------------------#
#	Esoh trimmomatic	    #
#-----------------------------------#

function usage() {
	printf "Usage: %s [P|N] [ options ]\n" $(basename $0);
	echo -e """
		Trim fastq files

                P - run parallel jobs
		N - no parallel jobs

		Options:
		-i,--idfile    <str>    :File containing SRA run accessions. [default: NULL]
		-p,--fqpath    <str>    :Path to fastq files. 
					:NB: Enter '.' for current directory 
					:[default: `pwd`/]
		-l,--leadx     <int>    :Number of bases to crop from the 5' end [default: 0]
		-t,--trailx    <int>    :Number of bases to crop from the 3' end [default: 0]
		-T,--threads   <int>    :Threads [default = 1]
		-h,--help               :Print this help message
	"""
}

if [ $# -lt 1 ]; then
    usage; 1>&2;
    exit 1;
fi

if [ $? != 0 ]; then
   echo "ERROR: Exiting..." 1>&2;
   exit 1;
fi

prog=`getopt -o "hi:p:l:t:T:" --long "help,idfile:,fqpath:,leadx:,trailx:,threads:" -- "$@"`

id=NULL
dname="`pwd`/"
leadx=0
trailx=0
t=1

eval set -- "$prog"

while true; do
    case "$1" in
      -i|--idfile) id="$2"; 
         if [[ "$2" == -* ]]; then 
            echo "ERROR: -i,--idfile must not begin with a '-'"; 1>&2; 
            exit 1; 
         fi; 
         shift 2 
         ;;
      -p|--fqpath) dname="$2/";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -p,--fqpath must not begin with a '-' "; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -l|--leadx) leadx="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -l,--leadx must not begin with a '-'"; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -t|--trailx) trailx="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -t,--trailx must not begin with a '-'"; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -T|--threads) t="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -T,--threads must not begin with a '-'"; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -h|--help) shift; usage; 1>&2; exit 1 ;;
      --) shift; break ;;
       *) shift; usage; 1>&2; exit 1 ;;
    esac
    continue
done

if [[ $id == NULL ]]; then
    echo "ERROR: -i,--idfile not provided! Terminating..."
    sleep 1;
    1>&2;
    exit 1;
else
echo -e """====================================
$(basename $0)		 kevin.esoh@students.jkuat.ac.ke (c) 2020

Option;		argument
idlist:		$id
fqpath:		$dname
leadx:		$leadx
trailx:		$trailx
threads:	$t
===================================="""
   which trimmomatic
   if [[ $? != 0 ]]; then
         echo "Could not locate trimmomatic"
         echo "Please install by one of the following means"
         echo "'sudo apt-get install trimmomatic'"
         echo "'conda install trimmomatic -c bioconda'"
         echo "Add bioconda channel if it doesn't already exist: 'conda config --add chennels bioconda'"
         1>&2;
         exit 1
   else
      mkdir -p paired unpaired
      awk -v d="${dname}" '{print d$1,d$2,"paired/"$1,"unpaired/"$1,"paired/"$1,"unpaired/"$1}' $id | sed 's/1.fastq.gz/fp.fastq.gz/2' |  sed 's/1.fastq.gz/fu.fastq.gz/2' | sed 's/1.fastq.gz/rp.fastq.gz/2' |  sed 's/1.fastq.gz/ru.fastq.gz/2' > trim.input.txt
      id=trim.input.txt
      if [[ "$1" == "P" ]]; then
         which parallel
         if [[ $? != 0 ]]; then
            echo "You have chosen to run jobs in parallel. However, GNU parallel does not appear to have been installed!"
            echo "Please install it and rerun your command"
            1>&2;
            exit 1;
         else
            echo "Your jobs will be run in parallel"
               cat $id | parallel --col-sep ' ' echo PE -phred33 {} ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:$leadx TRAILING:$trailx SLIDINGWINDOW:4:15 MINLEN:36 -threads $t | xargs -P10 -n15 trimmomatic
            rm $id
         fi
      else
         while read -r line; do
             trimmomatic PE -phred33 $line ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:$leadx TRAILING:$trailx SLIDINGWINDOW:4:15 MINLEN:36 -threads $t
         done < $id
         rm $id
      fi
   fi
fi
