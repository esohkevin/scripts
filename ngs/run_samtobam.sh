#!/bin/bash
cat samtobam.input.list | parallel --col-sep ' ' echo ./sam2bam.sh ../PlasmoDB-39_Pfalciparum3D7_Genome.fasta {1} {2} 10 | xargs -I {} -P 10 sh -c "{}"
