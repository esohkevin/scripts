#!/usr/bin/env bash

while getopt -l "invcf:,out:,nchr:,threads:,:help" -o "i:o:c:t:h" -- $@; do
   usage() { printf "Usage: %s: [-i|--invcf] [-o|--out] [-c|--nchr] [-t|--threads]\n" $(basename $0) 1>&2; exit 1; }
   case $@ in
     i|invcf) vcfin=$OPTARG;;
     o|out) oprfx=$OPTARG;;
     c|nchr) nchr=$OPTARG;;
     t|threads) thr=$OPTARG;;
     h|help) usage;;
     "") usage;;
     #*) usage;;
   esac
   echo -e "Arguments\nINPUT: $vcfin\nOUTPUT: $oprfx\nNUMCHR: $nchr\nTHREADS: $thr\n"
   sleep 1;
   echo "Starting job..."
   sleep 1;
   echo "<<< Progress... >>>"
   break
done
