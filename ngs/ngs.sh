#!/bin/bash

function usage() {
	printf "Usage: %s [ options ]\n" $(basename $0);
	echo -e """
		NGS Pipeline

		General Options:
		  -p,--fqpath    <str>    :Path to fastq files. 
					  :NB: Enter '.' for current directory 
		  			  :[default: `pwd`/]
                  -o,--out       <str>    :Output VCF file prefix [default: ngs ]
                  -T,--threads   <int>    :Threads [default: 1]
                  -h,--help               :Print this help message

		Trimmomatic Options:
		  -l,--leadx     <int>    :Number of bases to crop from the 5' end [default: 0]
		  -t,--trailx    <int>    :Number of bases to crop from the 3' end [default: 0]

		Alignment/Mapping Options:
		  -r,--ref       <str>	  :Reference FASTA sequence	
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

prog=`getopt -o "hp:o:l:t:r:T:" --long "help,fqpath:,out:,leadx:,trailx:,ref:,threads:" -- "$@"`

ref=NULL
dname="`pwd`/"
leadx=0
trailx=0
t=1
out="ngs"

eval set -- "$prog"

while true; do
    case "$1" in
      -p|--fqpath) dname="$2/";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -p,--fqpath must not begin with a '-' "; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -o|--out) out="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -o,--out must not begin with a '-' "; 1>&2;
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
      -r|--ref) ref="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -r,--ref must not begin with a '-'"; 1>&2;
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
if [[ "$ref" == NULL ]]; then
   echo "ERROR: -r,--ref not provided! Exiting..."; 1>&2;
   exit 1
fi

ls ${dname}*_1.fastq.gz > fwd.txt || ls ${dname}*_1.fq.gz > fwd.txt || ls ${dname}*_1.fastq > fwd.txt || ls ${dname}*_1.fq > fwd.txt
ls ${dname}*_2.fastq.gz > rev.txt || ls ${dname}*_2.fq.gz > rev.txt || ls ${dname}*_2.fastq > rev.txt || ls ${dname}*_2.fq > rev.txt
if [[ $? != 0 ]]; then
    echo "ERROR: No fastq file in the specified location Terminating..."
    sleep 1;
    1>&2;
    exit 1;
else
echo -e """========================================================
$(basename $0)		 kevin.esoh@students.jkuat.ac.ke (c) 2020

Option;		argument
fqpath:		$dname
ref:		$ref
leadx:		$leadx
trailx:		$trailx
threads:	$t
outFile:	$out
========================================================"""
    mkdir -p fastq paired unpaired aligned
    paste fwd.txt rev.txt | awk '{print $1,$2,$1}' | sed 's/_1.fastq.gz//2' > forward_reverse.txt
    awk '{print $1,$2}' forward_reverse.txt > fastq.input.txt
    id=forward_reverse.txt
    awk -v d="${dname}" '{print d$1,d$2,"paired/"$1,"unpaired/"$1,"paired/"$1,"unpaired/"$1}' $id | sed 's/1.fastq.gz/fp.fastq.gz/2' |  sed 's/1.fastq.gz/fu.fastq.gz/2' | sed 's/1.fastq.gz/rp.fastq.gz/2' |  sed 's/1.fastq.gz/ru.fastq.gz/2' > trim.input.txt
    awk '{print $3,$5}' trim.input.txt | sed 's/_1.fastq.gz/.sam/1' > align.input.txt
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
       id=fastq.input.txt; odr="fastq/"
       while read -r line; do
            echo "FastQC"
#           fastqc -t $t $line -o $odr
       done < $id
       
       id=trim.input.txt
       while read -r line; do
            echo "Trommomatic"
#           trimmomatic PE -phred33 $line ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:$leadx TRAILING:$trailx SLIDINGWINDOW:4:15 MINLEN:36 -threads $t
       done < $id
       
       bwa index $ref
       id=align.input.txt
       while read -r line; do
            echo "BWA"
#           bwa mem -t $t $ref $line -o aligned/$(awk '{print $1}' trim.input.txt | sed 's/_1.fastq.gz/.sam/1')
       done < $id
       for sam in aligned/$(awk '{print $1}' trim.input.txt | sed 's/_1.fastq.gz/.sam/1'); do
           samtools view -h ${sam} -O BAM -o ${sam/.sam/.bam}
           samtools sort -O BAM --reference $ref -@ $t -o ${sam/.sam/.sorted.bam} ${sam/.sam/.bam}
           echo ${sam/.sam/.sorted.bam}
       done > bam.list
       bcftools mpileup --min-MQ 2 --thread $t -f $ref -Oz -o out.vcf.gz -b bam.list
       bcftools index -f -t out.vcf.gz
       bcftools call -mv --threads $t -Oz -o ${out}.vcf.gz out.vcf.gz
       bcftools index -f -t ${out}.vcf.gz
       rm $id out.vcf.gz aligned/*.sam
    fi
fi
