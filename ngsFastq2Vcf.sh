#!/bin/bash

echo """
        #~@~# NGS Pipeline for Processing .fastq to .vcf #~@~#
              K. Esoh, kevin.esoh@students.jkuat.ac.ke
"""
######### Get recent versions of samtools and bcftools, and Trimmomatic ########
for tool in "samtools" "bcftools" "bwa" "bowtie2" "fastqc"; do    
    which $tool
	if [[ $? != "0"  ]]; then
	   if [[ $tool == "samtools" || $tool == "bcftools" ]]; then
	      echo "$tool is not installed!"
	      echo "If you have root access then install by running 'sudo apt-get install $tool' "
	      echo "If you don't have root access you can install a slightly lower version which can do most of the job (It may fail at some steps)"
	      read -p 'Would you like to install a lower version? [y|n] ' response
      	         if [[ $response == [Yy] ]]; then 
		    echo "Downloading... Please wait...";
		    wget https://sourceforge.net/projects/samtools/files/samtools/1.0/samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
		    bunzip2 samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
	            tar -xvf samtools-bcftools-htslib-1.0_x64-linux.tar
		    mkdir -p ~/bin
		    cp samtools-bcftools-htslib-1.0_x64-linux/bin/samtools ~/bin/
      		    samtools="samtools-bcftools-htslib-1.0_x64-linux/bin/samtools"
 		    bcftools="samtools-bcftools-htslib-1.0_x64-linux/bin/bcftools"
		 else 
		    echo "The pipeline can not work without $tool so it will terminate now. Thanks!"
	         fi
	    elif [[ $tool -eq "bwa" || $tool -eq "bowtie2" || $tool -eq "fastqc" ]]; then
	      echo "Please install $tool: 'sudo apt-get install $tool'"
	      echo "Or you visit the web to download"
      	    fi
	else
	   samtools="`which samtools`"
	   bcftools="`which bcftools`"
	fi
done

if [[ ! -f Trimmomatic-0.39/trimmomatic-0.39.jar ]]; then
	wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
	unzip Trimmomatic-0.39.zip
	chmod 755 Trimmomatic-0.39/trimmomatic-0.39.jar
else
	chmod 755 Trimmomatic-0.39/trimmomatic-0.39.jar
fi

cp -f Trimmomatic-0.39/adapters/* .
trimmomatic="Trimmomatic-0.39/trimmomatic-0.39.jar"

############################# Block-zip all fastq files #######################
if [[ -f "*.fq" ]]; then
   for fastqFile in ./*.fq; do
       mv ${fatsqFile} ${fatsqFile/.fq/.fastq}
   done
fi

if [[ -f "*.fastq" ]]; then
   for fastqFile in ./*.fastq ./*.fq; do
       bgzip ${fatsqFile}
   done
fi 

########################## Now initialize global parameters ####################
fastqBase=`for i in *1.fastq.gz *2.fastq.gz; do echo ${i/_*}; done | sort | uniq`
ref="PlasmoDB-39_Pfalciparum3D7_Genome.fasta"

########################### FastQc: Quality check #############################
#multiqc .
#fastqc
echo -e "\n#~@~# Trimmomatic #~@~#\n" 

for seqRead in ${fastqBase}; do
    java -jar $trimmomatic \
        PE \
        -phred33 \
        ${seqRead}_1.fastq.gz ${seqRead}_2.fastq.gz \
	${seqRead}_forward_paired.fastq.gz ${seqRead}_forward_unpaired.fastq.gz \
        ${seqRead}_reverse_paired.fastq.gz ${seqRead}_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
	-threads 10
done

###########################	   Make indexes     ############################
echo -e "\n#~@~# Indexing  Reference sequences #~@~#"
echo "` bwa index $ref `"
echo "` $samtools index -bc $ref ${ref/.fasta/.index} `"
echo "` $samtools faidx $ref `"

###########################      Alignment     #################################
echo -e "\n#~@~# Starting alignment #~@~#"

for readOrder in "forward_paired" "reverse_paired"; do
   for seqRead in ${fastqBase}; do
        
       bwa aln -t 50 $ref ${seqRead}_${readOrder}.fastq.gz > "${seqRead}"_"${readOrder}".sai
    
   done

done
echo -e "\nAll alignmnets completed: `date` "

#############	Combine paired-end reads to produce SAM files	 ###############
echo -e "\n#~@~# Mapping #~@~#"

for seqRead in ${fastqBase}; do
   bwa sampe \
	   	$ref ${seqRead}_forward_paired.sai \
		${seqRead}_reverse_paired.sai \
		${seqRead}_forward_paired.fastq.gz \
		${seqRead}_reverse_paired.fastq.gz \
		-f ${seqRead}_pe.sam
done
   echo -e "\nDone Mapping files: `date`"

########################## Convert SAM to BAM files ############################
echo -e "\n#~@~# Converting SAM to BAM files #~@~#"
echo "Please wait..."

for pairedReads in ${fastqBase}; do
    echo "`$samtools view -@ 10 -bS -T $ref ${pairedReads}_pe.sam -o ${pairedReads}.bam`"
done

################################ Sort BAM files ################################
echo -e "\n#~@~# Sorting BAM files #~@~#"
echo "Please wait..."

if [[ -f "bamlist.fofn" ]]; then
   rm bamlist.fofn
fi

for bam in ${fastqBase}; do
    echo ` $samtools sort -O bam -T temp.bam $bam.bam -o $bam.sorted.bam `
    echo $bam.sorted.bam >> bamlist.fofn
    echo `$samtools index ${bam}.sorted.bam`
done

########################### Merge sorted BAM files #############################
echo -e "\n#~@~# Merging Sorted BAM #~@~#"
echo "Please wait..."

$samtools merge -b bamlist.fofn Pf3D7.bam

###########################	 mpileup 	################################
echo -e "\n#~@~# Samtools mpileup - Variant Calling #~@~#\n"
echo "Please wait..."

$bcftools mpileup -d 100 --thread 20 -f $ref -Oz $i -o Pf3D7.vcf.gz -b bamlist.fofn

echo """
	Done Running all processes!: `date`
"""

