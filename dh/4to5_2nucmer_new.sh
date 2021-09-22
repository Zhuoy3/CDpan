#!/bin/bash
#SBATCH -J nucmer
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=200gb

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088 
do 
	cd /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge
	time nucmer -p /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs \
	/storage3/duh/pacbio/11_1/11_1.fasta \
	/storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs_rep_seq.fasta \
	-g 1000 -c 90 -l 500 -t 10
	show-coords -rcl \
	/storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.delta \
	> /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.delta.coords
done
