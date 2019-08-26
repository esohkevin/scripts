Conda pkgs
---
- conda build: Build pkgs


Terms
---
Conda channels: Sources repositories for conda pkgs. Conda-forge is one such channel. Also bioconda

Specifying a channel
```
conda install scipy --channel conda-forge
```

More than one channels
```
conda install numpy --channel conda-forge --channel bioconda
```

Suppress all default channels and only look from specified channel
```
conda install pandas --channel /home/user/path --overide-channels
```

Add alias to local channels
```
conda config --set channel_alias /home/esoh/bioTools
```

Install pkgs from local channel using alias name
```
conda install -c bioTools <pkg>
```


Bioinfo pkgs
---
**Install all packages without comment on availability by running the command-line below**
```
conda search <pkg>
conda install <pkg>=version [--channel | -c] <channel>  
```
Example
```
conda search eigensoft
conda install eigensoft=6.0.1 -c bioconda
conda install eigensoft=6.0.1 --channel bioconda
```

- Picard
- bcftools
- vcftools
- samtools
- bwa
- hisat2
- bowtie2
- plink1.9 (plink2 absent. Use linux-installed)
- emmax (Absent. Use linux-installed)
- bolt-lmm (Absent. Use linux-installed)
- eagle (v0.9.4.6 available. Use linux-installed v2.4)
- beagle (v4.1 available. Use)
- bedtools
- eigensoft
- fastqc
- multiqc
- tophat
- trimmomatic
- gatk3 (gatk4 not available. Use linux-installed)
> Download and unzip GATK4. Copy the gatk executable into home bin. Then run the line below to set gatk global environ variable
Preferably, place the line in _**~/.bashrc**_ file

```
export GATK_LOCAL_JAR=~/bioTools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar
```

- shapeit2
- qqman
