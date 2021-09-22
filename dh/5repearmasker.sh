#!/bin/bash
#SBATCH -J repeat
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 6

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088 
do
	cd /storage3/duh/pan_genome/Jiaxinghei/$i
	/home/duh/software/RepeatMasker/RepeatMasker -nolow -species pig /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.final.fa
	samtools faidx /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.final.fa
done
