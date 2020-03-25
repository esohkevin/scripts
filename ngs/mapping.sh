#!/bin/bash

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
   which bwa
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
      bwa index $ref
      if [[ "$1" == "P" ]]; then
         which parallel
         if [[ $? != 0 ]]; then
            echo "You have chosen to run jobs in parallel. However, GNU parallel does not appear to have been installed!"
            echo "Please install it and rerun your command"
            1>&2;
            exit 1;
         else
            echo "Your jobs will be run in parallel"
               cat $id | parallel --col-sep ' ' echo "mem -t $t $ref -R \\'@RG\\\tID:{3}\\\tPL:ILLUMINA\\\tSM:{3}\\' ${dname}{1} ${dname}{2} -o ${alnd}{3}.bwaPE.sam" | xargs -P10 -n10 bwa
            rm $id
         fi
      else
         while read -r line; do
             bwa mem -t $t $ref $line -o ${alnd}{3}.bwaPE.sam
         done < $id
         rm $id
      fi
   fi
fi


#export PATH=${PATH}:/BIODATA/StudentsHome/bin
if [[ $# == [34] ]]; then
    db="$HOME/db/"
    ref="$1"
    id="$2"
    dname="$3"
    t=$4
    mkdir -p aligned
    alnd="aligned/"
    bref="$(basename $ref)"
    paired_path="$2"
   
    conda activate cathyinformatics 
    
    echo -e "\n#~@~# Indexing  Reference sequence #~@~#"
    #hisat2-build $ref ${ref/.fa*/}
    #bowtie2-build $ref ${ref/.fa*/}
    bwa index $ref
    #novoindex $ref.nix $ref
    if [[ $# == 4 ]]; then
       n=$((50/$t))
       echo -e "\n#~@~# Starting alignment #~@~#"
       #cat $id | parallel --col-sep ' ' echo "-p $t -x ${ref} -1 $dname{1} -2 $dname{2} -S ${alnd}{3}.hisPE.sam" | xargs -P$n -n10 hisat2   
       #cat $id | parallel --col-sep ' ' echo "-p $t -x ${ref/.fa*/} -1 $dname{1} -2 $dname{2} -S ${alnd}{3}.bowPE.sam" | xargs -P$n -n10 bowtie2
       cat $id | parallel --col-sep ' ' echo "mem -t $t $ref -R \\'@RG\\\tID:{3}\\\tPL:ILLUMINA\\\tSM:{3}\\' ${dname}{1} ${dname}{2} -o ${alnd}{3}.bwaPE.sam" | xargs -P10 -n10 bwa
       for i in aligned/*.bwaPE.sam; do bgzip $i; done
       #cat $id | parallel echo "-d ${ref/.f*/.nix} -f $dname{/}_fp.fastq.gz $dname{/}_rp.fastq.gz -F ILMFQ -a -o SAM \\> ${alnd}{/}.novPE.sam" | xargs -P$n -n12 novoalign
    else
       #cat $id | parallel --col-sep ' ' echo "-x ${ref/.f*/} -1 $dname{1} -2 $dname{2} -S ${alnd}{3}.hisPE.sam" | xargs -P$n -n8 hisat2
       #cat $id | parallel --col-sep ' ' echo "-x ${ref/.*.bt2/} -1 $dname{1} -2 $dname{2} -S ${alnd}{3}.bowPE.sam" | xargs -P$n -n8 bowtie2
       cat $id | parallel --col-sep ' ' echo "mem -t 1 ${ref} -R \'@RG\\\tID:{3}\\\tPL:ILLUMINA\\\tSM:{3}\' ${dname}{1} ${dname}{2} -o ${alnd}{3}.bwaPE.sam" | xargs -P$n -n10 bwa
       #cat $id | parallel echo "-d ${ref/.f*/.nix} -f $dname{/}_fp.fastq.gz $dname{/}_rp.fastq.gz -F ILMFQ -a -o SAM \\> ${alnd}{/}.novPE.sam" | xargs -P$n -n12 novoalign

    fi
    #samtools index -bc $ref ${ref/.fasta/.index}

    echo -e "\nAll alignmnets completed: `date` "
    echo -e "Aligned SAM files can be found in $alnd"

else
    echo """
	Usage: ./mapping.sh <ref_seq> <idlist> <fpath> <threads>
	
	    ref: ref sequence (Provide full path)
	 idlist: File containing SRA run accessions. (NB: Must be in same path as fastq files)
	  fpath: Path to fastq files. (NB: Must end with a forward '/' slash)
	threads: Optional [default = 1]
	e.g.: File content may be of either format or both
	
	idlist Example: 3 columns, space delimited
	<forward> <reverse> <out-prefix>
	EEERROO1_fp.fastq.gz EEERR001_rp.fastq.gz EEERR001
    """
fi
