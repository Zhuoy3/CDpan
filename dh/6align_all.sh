#!/bin/bash
#SBATCH -J align
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 4

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088 
do
	mkdir -p /storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all
	cd /storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all
	bowtie2 \
	-x /storage3/duh/pan_genome/all/direct_merge/final/repeatmasker/index \
	-U /storage3/duh/pan_genome/Jiaxinghei/$i/unaligned_reads/singleUnmapped_R2.fq,/storage3/duh/pan_genome/Jiaxinghei/$i/unaligned_reads/singleUnmapped_R2.fq \
	-S /storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all/readContigAlignment.final.sam
done