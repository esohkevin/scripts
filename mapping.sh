#!/bin/bash

#export PATH=${PATH}:/BIODATA/StudentsHome/bin
if [[ $# == [34] ]]; then
    db="$HOME/db/"
    ref="$1"
    id="$2"
    dname="$3"
    t=$4
    mkdir -p ${dname}../aligned
    alnd="${dname}../aligned/"
    
    echo -e "\n#~@~# Indexing  Reference sequence #~@~#"
    #hisat2-build $ref ${ref/.fa*/}
    #bowtie2-build $ref ${ref/.fa*/}
    #bwa index $ref
    #novoindex $ref.nix $ref
    if [[ $# == 4 ]]; then
       n=$((50/$t))
       echo -e "\n#~@~# Starting alignment #~@~#"
       #cat $id | parallel echo "-p $t -x ${ref} -1 $dname{/}_fp.fastq.gz -2 $dname{/}_rp.fastq.gz -S ${alnd}{/}.hisPE.sam" | xargs -P$n -n10 hisat2   
       cat $id | parallel echo "-p $t -x ${ref/.fa*/} -1 $dname{/}_fp.fastq.gz -2 $dname{/}_rp.fastq.gz -S ${alnd}{/}.bowPE.sam" | xargs -P$n -n10 bowtie2
       cat $id | parallel echo "-t $t ${ref} $dname{/}_fp.fastq.gz $dname{/}_rp.fastq.gz -o ${alnd}{/}.bwaPE.sam" | xargs -P$n -n7 bwa mem
       #cat $id | parallel echo "-d ${ref/.f*/.nix} -f $dname{/}_fp.fastq.gz $dname{/}_rp.fastq.gz -F ILMFQ -a -o SAM \\> ${alnd}{/}.novPE.sam" | xargs -P$n -n12 novoalign
    else
       #cat $id | parallel echo "-x ${ref/.f*/} -1 $dname{/}_fp.fastq.gz -2 $dname{/}_rp.fastq.gz -S ${alnd}{/}.hisPE.sam" | xargs -P$n -n8 hisat2
       cat $id | parallel echo "-x ${ref/.*.bt2/} -1 $dname{/}_fp.fastq.gz -2 $dname{/}_rp.fastq.gz -S ${alnd}{/}.bowPE.sam" | xargs -P$n -n8 bowtie2
       cat $id | parallel echo "-t 1 ${ref} $dname{/}_fp.fastq.gz $dname{/}_rp.fastq.gz -o ${alnd}{/}.bwaPE.sam" | xargs -P$n -n7 bwa mem
       #cat $id | parallel echo "-d ${ref/.f*/.nix} -f $dname{/}_fp.fastq.gz $dname{/}_rp.fastq.gz -F ILMFQ -a -o SAM \\> ${alnd}{/}.novPE.sam" | xargs -P$n -n12 novoalign

    fi
    #samtools index -bc $ref ${ref/.fasta/.index}

    echo -e "\nAll alignmnets completed: `date` "

else
    echo """
	Usage: ./mapping.sh <ref_seq> <idlist> <fpath> <threads>
	
	    ref: ref sequence 
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
