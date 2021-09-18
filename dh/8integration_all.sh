#!/bin/bash
#SBATCH -J test
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 4

module load Perl
for breed in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S0352mu AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu; do
	##switch to destination file
	mkdir /storage3/duh/pan_genome/Anqing/$breed/link_new
	cd /storage3/duh/pan_genome/Anqing/$breed/link_new

	cp /storage3/duh/pan_genome/Anqing/$breed/bowtie2/all/mateLinks.txt ./

	mkdir ./1
	mkdir ./2
	mkdir ./3
	mkdir ./4

	##archive the contig name and length
	cp /storage3/duh/pan_genome/all/direct_merge/final/repeatmasker/dispensable_genome.fasta.masked.fai ./

	awk '{print $1,$2}' ./dispensable_genome.fasta.masked.fai \
		>./contig.name

	##拷贝操作脚本
	cp /storage3/duh/pan_genome/scripts/link_scripts/*.pl ./

	##first step(此步产生两个不同的文件，1文件为contig的两端都未有比对上参考基因组的，2为存在比对上参考基因组的)
	perl ./1do.pl ./contig.name ./mateLinks.txt

	##second step
	cat ./1/2to4.name | while read line; do
		contig=$(echo $line | awk '{print $1}')
		awk '{print $6}' ./4/$contig.link |
			sort - | uniq -c - >./4/$contig.chr
		sed -i 's/^[ ]*//' ./4/$contig.chr
		perl ./2do.pl ./4/$contig.chr $contig
	done 

	##third step
	mkdir -p ./2/2a
	mkdir -p ./2/2b/left
	mkdir -p ./2/2b/right

	cat ./3.contig.name | while read line; do
		ss=$(echo $line | awk '{print $1}')
		sort -u ./3/$ss.link >./3/$ss.link.new
	done

	awk 'NR==FNR{a[$1]=$2;next}{print $1,a[$1],$2}' ./contig.name ./3.contig.name >./3.contig.name.length

	cat ./3.contig.name.length |
		while read line; do
			contig=$(echo $line | awk '{print $1}')
			length=$(echo $line | awk '{print $2}')
			chr=$(echo $line | awk '{print $3}')
			perl ./3do.pl -n $contig \
				-i ./3/$contig.link.new -l $length -c $chr
		done

	##fourth step
	for i in $(ls ./4/); do
		if [ -s ./4/$i ]; then
			contig=$(echo $i | awk -F'.' '{print $1}')
		fi
		echo $contig >>./4.contig.name
	done

	awk 'NR==FNR{a[$1]=$2;next}{print $1,a[$1],$2}' ./contig.name ./4.contig.name >./4.contig.name.length

	cat ./4.contig.name.length |
		while read line; do
			contig=$(echo $line | awk '{print $1}')
			length=$(echo $line | awk '{print $2}')
			perl ./4do.pl -n $contig \
				-i ./4/$contig.link -l $length
		done

	##fifth step
	for i in $(ls ./4/); do
		if [ -s ./4/$i ]; then
			contig=$(echo $i | awk -F'.' '{print $1}')
		fi
		echo $contig >>./1/4.name
	done

	for i in $(ls ./3/); do
		if [ -s ./3/$i ]; then
			contig=$(echo $i | awk -F'.' '{print $1}')
		fi
		echo $contig >>./1/3.name
	done

	rm -rf ./1/2to4.name

	awk 'NR==FNR{a[$1]=$2;next}{print $1,a[$1],$2}' ./contig.name ./1/4.name >./1/4.name.re
	awk 'NR==FNR{a[$1]=$2;next}{print $1,a[$1],$2}' ./contig.name ./1/3.name >./1/3.name.re

	mv ./1/4.name.re ./1/4.name
	mv ./1/3.name.re ./1/3.name
done
