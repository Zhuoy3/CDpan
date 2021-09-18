#!/bin/bash
###
 # @Description:
 # @Author: Zhuo Yue
 # @Date: 2021-07-15 23:00:43
 # @LastEditors: Zhuo Yue
 # @LastEditTime: 2021-08-24 11:35:55
 # @Calls:
 # @Called By:
 # @FilePath: \CDpan\dh\4judge.sh
###
#PBS -N Anqing
#PBS -q batch

for i in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu AQF-S0352mu
do
	cd /diskd/duh/domestic/pan_genome/Anqing/$i
	perl /diskd/duh/domestic/pan_genome/extract_id.pl \
	/diskd/duh/domestic/pan_genome/Anqing/$i/centrifuge.output \
	/diskd/duh/domestic/pan_genome/centrifuge/Chordata.taxid

	perl /diskd/duh/domestic/pan_genome/extract_sequence.pl \
	/diskd/duh/domestic/pan_genome/Anqing/$i/keep_id.list \
	/diskd/duh/domestic/pan_genome/Anqing/$i/large_1000.fasta
done
