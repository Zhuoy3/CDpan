#!/bin/bash
#SBATCH -J mmseqs
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 10

source /home/duh/miniconda3/bin/activate mmseqs2
for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088 
do 
	cd /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge
	time mmseqs easy-linclust \
	/storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.fa \
	/storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs \
	/storage3/duh/tmp --threads 10 --cov-mode 1 -c 0.9 --min-seq-id 0.9 --cluster-mode 2
done
