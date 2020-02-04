#!/usr/bin/env bash

while getopt -l "vcf:,oprefix:nchr:,threads::help" -o "v:o:c:t::h" -- $@; do
   usage() { printf "\nUsage: %s: [-v|--vcf] [-o|--oprefix] [-c|--nchr] [-t|--threads]\n" $(basename $0) 1>&2; exit 1; }
   case $@ in
     v) vcfin=$OPTARG;;
     o) oprfx=$OPTARG;;
     c) nchr=$OPTARG;;
     t) thr=$OPTARG;;
     h) usage;;
     "") usage;;
     #*) usage;;
   esac
done
