#!/bin/bash

zcat chr_updated.vcf.gz | \
	sed 's/Pf3D7_01_v3/1/g' | \
	sed 's/Pf3D7_02_v3/2/g' | \
	sed 's/Pf3D7_03_v3/3/g' | \
	sed 's/Pf3D7_04_v3/4/g' | \
	sed 's/Pf3D7_05_v3/5/g' | \
	sed 's/Pf3D7_06_v3/6/g' | \
	sed 's/Pf3D7_07_v3/7/g' | \
	sed 's/Pf3D7_08_v3/8/g' | \
	sed 's/Pf3D7_09_v3/9/g' | \
	sed 's/Pf3D7_10_v3/10/g' | \
	sed 's/Pf3D7_11_v3/11/g' | \
	sed 's/Pf3D7_12_v3/12/g' | \
	sed 's/Pf3D7_13_v3/13/g' | \
	sed 's/Pf3D7_14_v3/14/g'  > chrIdUpdated.vcf

bgzip chrIdUpdated.vcf
tabix -p vcf chrIdUpdated.vcf.gz
