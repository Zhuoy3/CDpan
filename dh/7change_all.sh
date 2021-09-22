#!/bin/bash
#SBATCH -J test
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 4

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088 
do
	samtools view -h -F 256 \
	/storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all/readContigAlignment.final.sam |\
	samtools sort - -n -O bam |bedtools bamtobed -i stdin |\
	awk '{OFS="\t"}{print $4,$1,$6,$2,$3}' |\
	sort > /storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all/readContigAlignment.txt

	samtools view -H \
	/storage3/duh/pan_genome/Jiaxinghei/$i/unaligned_reads/Sus_11_1_Links.bam |\
	cat - <(awk 'FNR==NR{main[$1]=$0;next} $1 in main {print main[$1]}' <(samtools view /storage3/duh/pan_genome/Jiaxinghei/$i/unaligned_reads/Sus_11_1_Links.bam) /storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all/readContigAlignment.txt) |\
	samtools sort -n -O bam |bedtools bamtobed -i stdin |\
	awk '{OFS="\t"}{print $4,$1,$6,$2,$3}' |sed -e 's/\/[1-2]//g' |\
	sort > /storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all/matchedMates.txt

	join -j 1 \
	/storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all/readContigAlignment.txt \
	/storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all/matchedMates.txt \
	> /storage3/duh/pan_genome/Jiaxinghei/$i/bowtie2/all/mateLinks.txt
done
