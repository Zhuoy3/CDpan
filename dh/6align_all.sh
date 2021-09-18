#!/bin/bash
###
# @Description:
# @Author: Zhuo Yue
# @Date: 2021-07-15 23:00:43
# @LastEditors: Zhuo Yue
# @LastEditTime: 2021-08-30 11:13:00
# @Calls:
# @Called By:
# @FilePath: \CDpan\dh\6align_all.sh
###
#SBATCH -J align
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 6

for i in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S0352mu AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu; do
	mkdir -p /storage3/duh/pan_genome/Anqing/$i/bowtie2/all
	cd /storage3/duh/pan_genome/Anqing/$i/bowtie2/all
	bowtie2 \
		-x /storage3/duh/pan_genome/all/direct_merge/final/repeatmasker/index \
		-U /storage3/duh/pan_genome/Anqing/$i/unaligned_reads/singleUnmapped_R2.fq,/storage3/duh/pan_genome/Anqing/$i/unaligned_reads/singleUnmapped_R2.fq \
		-S /storage3/duh/pan_genome/Anqing/$i/bowtie2/all/readContigAlignment.final.sam -p 6
done
