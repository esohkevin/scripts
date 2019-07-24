#!/bin/bash

if [[ $# != 3 ]]; then
   echo """
	Usage: ./runScan.sh <.hap+.map-root> <chr#> <pop-name>
"""
else
   Rscript runScan.R $1 $2 $3

   echo """
        Done Scanning Genome! Results saved in iHSResult.txt, frqResults.txt, and chr11Signals.txt
        Please check the 'chr11Signals.txt' file for genomic positions of interest and run the 'runEHH.sh' script for EHH
"""

fi
