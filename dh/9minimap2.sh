#!/bin/bash
#SBATCH -J minimap2
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 2

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088 
do
	mkdir /storage3/duh/pan_genome/Jiaxinghei/$i/minimap2
	cd /storage3/duh/pan_genome/Jiaxinghei/$i/minimap2
	/home/duh/software/minimap2/minimap2 -x asm10 \
	/storage3/duh/pan_genome/all/direct_merge/final/dispensable_genome.fasta \
	/storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.final.fa \
	> /storage3/duh/pan_genome/Jiaxinghei/$i/minimap2/aln.paf
done
