#!/bin/bash
#SBATCH -J nucmer
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=200gb

for i in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S0352mu AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu
do
	cd /storage3/duh/pan_genome/Anqing/$i/centrifuge
	time nucmer -p /storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs \
	/storage3/duh/pacbio/11_1/11_1.fasta \
	/storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs_rep_seq.fasta \
	-g 1000 -c 90 -l 500 -t 30
	show-coords -rcl \
	/storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.delta \
	> /storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.delta.coords
done
