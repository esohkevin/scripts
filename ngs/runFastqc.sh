#!/bin/bash

#!/bin/bash

function usage() {
  printf "Usage: %s [p/n] [ Options ]\n" $(basename $0)
  echo -e """
	  FastQC 
	 
	  NB: The first argument must be either 'p' or 'n' to specify parallel jobs or not
	  ________________________________________________________________________________________
	   p					:Run parallel jobs [Requires GNU parallel]
	   n					:Run jobs in serial (No parallele jobs)
	  ________________________________________________________________________________________

        Options:
	  -i,--idfile	<str>	:File containing fastq file names, two columns, where column 1
				 is forward run files and column 2 is the reverse run files.

				 Examples: (NB: space or tab delimited)
					ERRAOOOOO1_R1.fastq.gz ERRAOOOOO1_R2.fastq.gz
                                        ERRAOOOOO2_R1.fastq.gz ERRAOOOOO2_R2.fastq.gz

	  -f,--fqpath   <str>	:Path to fastq files [default: `pwd`/]
	  -t,--threads  <int>	:Number of threads [default: 1]
  """
}

if [ $# -lt 1 ]; then
   usage; 1>&2;
   exit 1;
fi

if [ $? != 0 ]; then
   echo "ERROR: Exiting..."; 1>&2;
   exit 1;
fi

prog=`getopt -o "hi:f:t:" --long "help,idfile:,fqpath:,threads:" -- "$@"`

id=NULL
dname="`pwd`/"
thr=1

eval set -- "$prog"

while true; do
   case "$1" in
      -i|--idfile) id="$2"
         if [[ "$2" == -* ]]; then
            echo "ERROR: -i,--idfile must not begin with a '-'"; 1>&2;
            exit 1;
         fi
         shift 2
         ;;
      -f|--fqpath) dname="$2/"
         if [[ "$2" == -* ]]; then
            echo "ERROR: -f,--fqpath must not begin with a '-'"; 1>&2;
            exit 1;
         fi
         shift 2;
         ;;
      -t|--threads) thr="$2"
         if [[ "$2" == -* ]]; then
            echo "ERROR: -t,--threads must not start with a '-'"; 1>&2;
            exit 1;
         fi
         shift 2;
         ;;
      -h|--help) shift; usage; 1>&2; exit 1 ;;
      --) shift; break ;;
       *) shift; usage; 1>&2; exit 1 ;;
   esac
   continue
done

if [[ "$id" == NULL ]]; then
   echo "ERROR: -i,--idfile not provided! Terminating..."; 1>&2;
   sleep 1;
   exit 1;
else
echo -e """==============================================
$(basename $0)		kevin.esoh@students.jkuat.ac.ke (c) 2020
INPUT DETAILS
Option:			argument
idfile:			$id
fqpath:			$dname
threads:		$thr
=============================================="""
mkdir -p fastq
odr="fastq/"
awk '{print $1,$2}' project.ids > forward.reverse.ids
id=forward.reverse.ids
   if [[ "$1" == "p" ]]; then
      which parallel
      if [[ $? != 0 ]]; then
         echo "You have chosen to run jobs in parallel. However, GNU parallel does not appear to have been installed!"
         echo "Please install it and rerun your command"
         1>&2;
         exit 1;
      else
         echo "Your jobs will be run in parallel"
            cat $id | parallel --col-sep ' ' echo "-t $thr $dname{1} $dname{2} -o $odr" | xargs -P5 -n6 fastqc
         rm $id
      fi
   else
      while read -r line; do
         fastqc -t $thr $line -o $odr
      done < $id
   fi
fi
