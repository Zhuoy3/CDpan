#!/bin/bash
#SBATCH -J test
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 4

for i in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S0352mu AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu
do
	samtools view -h -F 256 \
	/storage3/duh/pan_genome/Anqing/$i/bowtie2/all/readContigAlignment.final.sam |\
	samtools sort - -n -O bam |bedtools bamtobed -i stdin |\
	awk '{OFS="\t"}{print $4,$1,$6,$2,$3}' |\
	sort > /storage3/duh/pan_genome/Anqing/$i/bowtie2/all/readContigAlignment.txt

	samtools view -H \
	/storage3/duh/pan_genome/Anqing/$i/unaligned_reads/Sus_11_1_Links.bam |\
	cat - <(awk 'FNR==NR{main[$1]=$0;next} $1 in main {print main[$1]}' <(samtools view /storage3/duh/pan_genome/Anqing/$i/unaligned_reads/Sus_11_1_Links.bam) /storage3/duh/pan_genome/Anqing/$i/bowtie2/all/readContigAlignment.txt) |\
	samtools sort -n -O bam |bedtools bamtobed -i stdin |\
	awk '{OFS="\t"}{print $4,$1,$6,$2,$3}' |sed -e 's/\/[1-2]//g' |\
	sort > /storage3/duh/pan_genome/Anqing/$i/bowtie2/all/matchedMates.txt

	join -j 1 \
	/storage3/duh/pan_genome/Anqing/$i/bowtie2/all/readContigAlignment.txt \
	/storage3/duh/pan_genome/Anqing/$i/bowtie2/all/matchedMates.txt \
	> /storage3/duh/pan_genome/Anqing/$i/bowtie2/all/mateLinks.txt
done
