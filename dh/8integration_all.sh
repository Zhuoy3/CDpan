#!/bin/bash
#SBATCH -J test
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 4

module load Perl
for breed in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088; do
	##switch to destination file
	mkdir /storage3/duh/pan_genome/Jiaxinghei/$breed/link_new
	cd /storage3/duh/pan_genome/Jiaxinghei/$breed/link_new

	cp /storage3/duh/pan_genome/Jiaxinghei/$breed/bowtie2/all/mateLinks.txt ./

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
			echo $contig >>./4.contig.name
		else
			unset contig
		fi
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
			echo $contig >>./1/4.name
		else
			unset contig
		fi
	done

	for i in $(ls ./3/); do
		if [ -s ./3/$i ]; then
			contig=$(echo $i | awk -F'.' '{print $1}')
			echo $contig >>./1/3.name
		else
			unset contig
		fi
	done

	cat ./1/3.name | while read line; do
		contig=$(echo $line | awk '{print $1}')
		awk '{print $6}' ./3/$contig.link |
			sort - | uniq -c - | sed 's/^[ ]*//' - | sort -nrk 1 - | head -n 1 - >./3/$contig.chr
		chrpart3=$(awk '{print $2}' ./3/$contig.chr)
		printf "$line $chrpart3\n" >>./1/3.name.re
	done

	rm -rf ./1/2to4.name

	awk 'NR==FNR{a[$1]=$2;next}{print $1,a[$1],$2}' ./contig.name ./1/4.name >./1/4.name.re

	mv ./1/4.name.re ./1/4.name
	mv ./1/3.name.re ./1/3.name

	for i in $(ls ./4/); do
		if [ ! -s ./4/$i ]; then
			contig=$(echo $i | awk -F'.' '{print $1}')
			echo $contig >>./1/5.name
		else
			unset contig
		fi
	done

	rm -rf ./4
	mkdir ./4
done
