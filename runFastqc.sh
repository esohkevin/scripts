#!/bin/bash

if [[ $# == [23] ]]; then
   id="$1"
   dname="$2"
   t=$3
   mkdir -p ${dname}fastqced
   dr="${dname}fastqced/"
   if [[ $# == 3 ]]; then
      n="$((50/$t))" 
      cat $id | parallel echo "-t $t $dname{/}_1.fastq.gz $dname{/}_2.fastq.gz -o $dr" | xargs -P$n -n6 fastqc
   else
      cat $id | parallel echo "-t 1 {/}_1.fastq.gz {/}_2.fastq.gz -o $dr" | xargs -P5 -n6 fastqc
   fi
   echo "Done! All results saved in '$dr'"
else
   echo -e """
	Usage: runFastqc <idlist> <fpath> <threads>

	 idlist: File containing SRA run accessions. (NB: Must be in same path as fastq files)
          fpath: Path to fastq files. (NB: Must end with a forward '/' slash)
        threads: Optional [default = 1]

        e.g.: File content may be of either format or both
        
        ftp.sra.ebi.ac.uk/vol1/err/ERR182/007/ERR1823587
        ftp.sra.ebi.ac.uk/vol1/err/ERR182/008/ERR1823588
        ftp.sra.ebi.ac.uk/vol1/err/ERR182/009/ERR1823589
        ERR1823590
        ERR1823591
        ERR1823592
	
    """
fi
