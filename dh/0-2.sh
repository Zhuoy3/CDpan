#!/bin/bash
#PBS -N anqing
#PBS -q lightning
#PBS -l nodes=1:ppn=12

for i in AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu AQF-S0352mu
do
	mkdir /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i
	samtools fastq -@ 12 -f 12 /storage3/DATA_storage3/Anqingfeidong/$i/$i.recal.bam -1 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/mateUnmapped_R1.fq -2 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/mateUnmapped_R2.fq
	samtools fastq -@ 12 -f 68 -F 8 /storage3/DATA_storage3/Anqingfeidong/$i/$i.recal.bam > /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/R1_mateMapped.fq
	samtools fastq -@ 12 -f 132 -F 8 /storage3/DATA_storage3/Anqingfeidong/$i/$i.recal.bam > /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/R2_mateMapped.fq
	samtools view -@ 12 -f 8 -F 4 /storage3/DATA_storage3/Anqingfeidong/$i/$i.recal.bam > /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/Sus_11_1_Links.bam
done