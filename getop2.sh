#!/bin/bash

#set opt=""

function usage() {
	printf "Usage: %s [-i|--invcf] [-o|--out] [-n|--nchr] ...\n" $(basename $0)
	echo """
		This program takes in an parses arguments

		-i <str>	:VCF input file. Specify the path [Required]
		-o <str>	:Output file [default: invcf_base.out]
		-n <int>	:Number chromosomes you wish to run [default: 1]
		-t <int>	:Number of extra threads you wish to use [default: 1]
		-T <int>	:The iHS threshold you wish to set [default: computed on the fly]
		-h <str>	:Get this help message
	"""
}

while getopts "hi:o:n:t:T:" opt; do
    case "$opt" in
       h|':'|'?'|'*') usage; 1>&2; exit 1; ;;
       i) inv=${OPTARG} ;;
       o) out=${OPTARG} ;;
       n) nc=${OPTARG} ;;
       t) thr=${OPTARG} ;;
       T) thresh=${OPTARG} ;;
    esac
done

#while getopt -a "hi:o:n:t:T:" -- opt; do
#    case "$opt" in
#       h|help) usage; 1>&2; exit 1; ;;
#       i|invcf) inv=${OPTARG} ;;
#       o|out) out=${OPTARG} ;;
#       n|nchr) nc=${OPTARG} ;;
#       t|threads) thr=${OPTARG} ;;
#       T|thresh) thesh=${OPTARG} ;;
#    esac
#    shift $(( $OPTIND-1 ))
#done

sleep 1
echo "Arguments received"
echo -e "INPUT: $inv\nOUTPUT: $out\nNUMCHR: $nc\nTHREADS: $thr\nTHRESH: $thresh\n"
sleep 1
echo "Now preparing to run jobs..."
sleep 1
echo -e "Executing\e[32;5;5m...\e[0m"


