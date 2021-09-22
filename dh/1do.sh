#!/bin/bash
#PBS -N Jiaxinghei
#PBS -q lightning
#PBS -l nodes=1:ppn=12

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088
do
	mkdir /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i
	samtools fastq -@ 12 -f 12 /storage2/liujf/Jiaxinghei/$i/$i.recal.bam -1 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/mateUnmapped_R1.fq -2 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/mateUnmapped_R2.fq
	samtools fastq -@ 12 -f 68 -F 8 /storage2/liujf/Jiaxinghei/$i/$i.recal.bam > /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/R1_mateMapped.fq
	samtools fastq -@ 12 -f 132 -F 8 /storage2/liujf/Jiaxinghei/$i/$i.recal.bam > /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/R2_mateMapped.fq
	samtools view -@ 12 -f 8 -F 4 /storage2/liujf/Jiaxinghei/$i/$i.recal.bam > /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/Sus_11_1_Links.bam
done