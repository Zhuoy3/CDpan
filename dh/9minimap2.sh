#!/bin/bash
#SBATCH -J minimap2
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 2

for i in BM01 BM02 BM03 BM06 BM09 BM19 BM25 BM27 BM29 BM30 BM33 BM35 BM38 BM40 BM41 BM45 BM47 BM48 BM52 BM55 BM57 BM58
do
	mkdir /storage3/duh/pan_genome/Bama/$i/minimap2
	cd /storage3/duh/pan_genome/Bama/$i/minimap2
	/home/duh/software/minimap2/minimap2 -x asm10 \
	/storage3/duh/pan_genome/all/direct_merge/final/dispensable_genome.fasta \
	/storage3/duh/pan_genome/Bama/$i/centrifuge/filtered.mmseqs.final.fa \
	> /storage3/duh/pan_genome/Bama/$i/minimap2/aln.paf
done
