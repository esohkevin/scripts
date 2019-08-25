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
- Picard
- bcftools
- vcftools
- samtools
- bwa
- hisat2
- bowtie2
- plink1.9/2
- emmax
- bolt-lmm
- eagle
- beagle
- gatk
- bedtools
- eigensoft
- fastqc
- multiqc
- tophat
