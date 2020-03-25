#!/usr/bin/env bash

if [[ $# != 4 ]]; then
   echo -e "\nUsage: ./sam2bam.sh <reference> <insam> <outbam> <threads>\n"
else
   ref=$1
   insam=$2
   outbam=$3
   thr=$4
   samtools view -@ $thr -h -q 2 -F 4 -bT $ref $insam | samtools sort - > $outbam
fi
