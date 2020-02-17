Random scripts
---
![samtools filtering options](https://github.com/esohkevin/scripts/blob/master/samFilter.png)

- Python
- Shell
- R
- etc

For generating ssh keys, click [here](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/)

Conda pkgs
---
- conda build: Build pkgs


Terms
---
Conda channels: Source repositories for conda pkgs. Conda-forge is one such channel. Also bioconda

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
conda install pandas --channel file:///home/user/path --overide-channels
```

Add alias to local channels
```
conda config --set channel_alias file:///home/esoh/bioTools
```

Install pkgs from local channel using alias name
```
conda install -c bioTools <pkg>
```

Make a copy of all installed packages in current env and use to create a new identical env
```
conda list --explicit > spec-file.txt
conda create --name newenv --file spec-file.txt
```

Create identical env from an env.yml file
```
conda env create [-f|--file] environment.yml 
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
- admixture
- finestructure
- hierfstat (R)
- pegas (read VCF in R)
- adegenet (R)
- gcta (NA)
- snptest
- hapFlk (selection scan)
- phangorn (genetic trees)
- vtools (bioconda) ***NGS, Exome Seq Analysis***
- modeller

```
conda config --add channels salilab
KEY_MODELLER='XXX'
conda install modeller
```
- Pymol

```
conda config --prepend channels schrodinger
conda install -c schrodinger pymol
```
