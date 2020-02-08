#!/bin/bash

#-------------------------------------- GETOPTS --------------------------------------
#function usage() {
#	printf "Usage: %s [-i invcf] [-o outfile] [-n #chr] ...\n" $(basename $0)
#	echo """
#		This program takes in an parses arguments
#
#		-i,-invcf	<str>		:VCF input file. Specify the path [Required]
#		-o,-out 	<str>		:Output file [default: invcf_base.out]
#		-n,-nchr 	<int>		:Number chromosomes you wish to run [default: 1]
#		-t,-threads 	<int>		:Number of extra threads you wish to use [default: 1]
#		-T,-thresh 	<int>		:The iHS threshold you wish to set [default: computed on the fly]
#		-h,-help 	<str>		:Show this help message
#	"""
#}

#while getopts "hi:o:n:t:T:" opt; do
#    case "$opt" in
#       h|help) usage; 1>&2; exit 1 ;;
#       i|invcf) inv=${OPTARG} ;;
#       o|out) out=${OPTARG} ;;
#       n|nchr) nc=${OPTARG} ;;
#       t|threads) thr=${OPTARG} ;;
#       T|thresh) thresh=${OPTARG} ;;
#       '?') usage; 1>&2; exit 1 ;;
#       '*') usage; 1>&2; exit 1 ;;
#    esac
#done

#-------------------------------------- GETOPT ---------------------------------------
function usage() {
        printf "Usage: %s [-i invcf] [-o outfile] [-n #chr] ...\n" $(basename $0)
        echo """
                This program takes in an parses arguments

                -i,--invcf       <str>           :VCF input file. Specify the path [Required]
                -o,--out         <str>           :Output file [default: invcf_base.out]
                -n,--nchr        <int>           :Number chromosomes you wish to run [default: 1]
                -t,--threads     <int>           :Number of extra threads you wish to use [default: 1]
                -T,--thresh      <int>           :The iHS threshold you wish to set [default: computed on the fly]
                -h,--help        <NULL>          :Show this help message
		-v,--verbose     <NULL>		 :Verbose
        """
}

temp=`getopt -o "hi:o:n:t:T:v" -l "help,invcf:,out:,nchr:,threads:,thresh:,verbose" -- "$@"`

if [ $? != 0 ]; then 
   echo "Terminating..." 1>&2; 
   exit 1; 
fi

VERBOSE=false
DEBUG=false
eval set -- "$temp"

while true; do
    case "$1" in
       -i|--invcf) inv="$2"; shift 2;;
       -o|--out) out="$2"; shift 2;;
       -n|--nchr) nc="$2"; shift 2;;
       -t|--threads) thr="$2"; shift 2;;
       -T|--thresh) thresh="$2"; shift 2 ;;
       -v|--verbose) VERBOSE=true; shift ;;
       -h|--help) shift; usage; 1>&2; exit 1 ;;
       --) shift; break ;;
       *) shift; usage; 1>&2; exit 1 ;;
    esac
done

#-------------------------------------- PRINT ARGS -----------------------------------
sleep 1
echo "Arguments received"
echo -e "INPUT: $inv\nOUTPUT: $out\nNUMCHR: $nc\nTHREADS: $thr\nTHRESH: $thresh\n"
sleep 1
echo "Now preparing to run jobs..."
sleep 1
echo -e "Executing\e[32;5;5m...\e[0m"

