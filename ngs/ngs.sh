#!/bin/bash

function usage() {
	printf "Usage: %s [ fq|trim|map ] [ options ]\n" $(basename $0);
	echo -e """
		NGS Pipeline - Paired End

		Enter 'fq' OR 'trim' OR 'map' OR all to run either FastQC OR Trimmomatic OR alignment/mapping or all steps.
		Enter any two to run only those two steps

		General Options:
		  -p,--fqpath    <str>    :Path to fastq files. NB: Make sure all the files are paired i.e. forward and reverse 
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

for i in $( (ls ${dname}*_1.fastq.gz || ls ${dname}*_1.fq.gz || ls ${dname}*_1.fastq || ls ${dname}*_1.fq) && (ls ${dname}*_R1*.fastq.gz || ls ${dname}*_R1*.fq.gz || ls ${dname}*_R1*.fastq || ls ${dname}*_R1*.fq) ); do
    if [[ -f "$i" ]]; then
        basename $i;
    else
	echo "ERROR: No fastq file file in the specified location Terminating..."
	1>&2;
	exit 1
    fi
done > fwd.txt

for i in $( (ls ${dname}*_2.fastq.gz || ls ${dname}*_2.fq.gz || ls ${dname}*_2.fastq || ls ${dname}*_2.fq) && (ls ${dname}*_R2*.fastq.gz || ls ${dname}*_R2*.fq.gz || ls ${dname}*_R2*.fastq || ls ${dname}*_R2*.fq) ); do
    if [[ -f "$i" ]]; then
       basename $i;
    else
        echo "ERROR: No fastq file found in the specified location Terminating..."
        1>&2;
        exit 1
    fi
done > rev.txt

if [[ $? != 0 ]]; then
    echo "ERROR: No fastq file in the specified location Terminating..."
    sleep 1;
    1>&2;
    exit 1;
else
    echo -e """
               ===================================================================
               NGS Pipeline		 kevin.esoh@students.jkuat.ac.ke (c) 2020
               -------------------------------------------------------------------
               Option;		argument
               fqpath:		$dname
               reference:	$ref
               leading:		$leadx
               trailing:	$trailx
               threads:		$t
               outFile:		${out}.vcf.gz
               ===================================================================
               Starting NGS Pipeline. Please wait...
    """
    mkdir -p fastq paired unpaired aligned
    paste fwd.txt rev.txt | awk '{print $1,$2}' > forward_reverse.txt
    for i in $(cat fwd.txt); do if [[ ! -f "paired/${i}_fp.fq.gz" ]]; then cp ${dname}${i} paired/${i}_fp.fq.gz; fi; done
    for i in $(cat rev.txt); do if [[ ! -f "paired/${i}_rp.fq.gz" ]]; then cp ${dname}${i} paired/${i}_rp.fq.gz; fi; done
    awk -v d="${dname}" '{print d$1,d$2}' forward_reverse.txt > fastq.input.txt
    awk -v d="${dname}" '{print d$1,d$2,"paired/"$1"_fp.fq.gz","unpaired/"$1"_fu.fq.gz","paired/"$2"_rp.fq.gz","unpaired/"$2"_ru.fq.gz"}' forward_reverse.txt > trim.input.txt
    awk '{print $3,$5,$3}' trim.input.txt | sed 's/_fp.fq.gz/.sam/2' | sed 's/paired\///3' |  awk '{print $1,$2,"-o","aligned/"$3}' > align.input.txt
    rm fwd.txt rev.txt
    while true; do
      case "$1" in
         fq)
           id=fastq.input.txt; odr="fastq/"
           while read -r line; do
               echo "FastQC"
               fastqc -t $t $line -o $odr
           done < $id
           shift
           ;;
         trim)
           id=trim.input.txt
           while read -r line; do
               echo "Trommomatic"
               trimmomatic PE -phred33 $line ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:$leadx TRAILING:$trailx SLIDINGWINDOW:4:15 MINLEN:36 -threads $t
           done < $id
           shift
           ;;
         map)
           if [[ "$ref" == NULL ]]; then
              echo "ERROR: -r,--ref not provided! Exiting..."; 1>&2;
              exit 1
	   elif [[ ! -f "${ref}.bwt" ]]; then
		bwa index $ref
           fi
           id=align.input.txt
           while read -r line; do
               echo "BWA"
               bwa mem -t $t $ref $line
           done < $id
           for sam in $(awk '{print $4}' align.input.txt); do
               samtools view -h ${sam} -O BAM -o ${sam/.sam/.bam}
               samtools sort -O BAM --reference $ref -@ $t -o ${sam/.sam/.sorted.bam} ${sam/.sam/.bam}
               echo ${sam/.sam/.sorted.bam}
               rm ${sam/.sam/.bam}
           done > bam.list
           bcftools mpileup --min-MQ 2 --thread $t -f $ref -Oz -o out.vcf.gz -b bam.list
           bcftools index -f -t out.vcf.gz
           bcftools call -mv --threads $t -Oz -o ${out}.vcf.gz out.vcf.gz
           bcftools index -f -t ${out}.vcf.gz
           rm $id out.vcf.gz aligned/*.sam
	   shift
           ;;
	 *) shift; 1>&2; exit 1 ;;
      esac
      continue
    done
fi
