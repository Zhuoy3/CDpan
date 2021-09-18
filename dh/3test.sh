#!/bin/bash
###
 # @Description:
 # @Author: Zhuo Yue
 # @Date: 2021-07-15 23:00:43
 # @LastEditors: Zhuo Yue
 # @LastEditTime: 2021-07-23 09:53:36
 # @Calls:
 # @Called By:
 # @FilePath: \CDpan\dh\3test.sh
###
#PBS -N test
#PBs -q batch

for i in AQ-2 AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu AQF-S0352mu
do
	cd /diskd/duh/domestic/pan_genome/Anqing/$i
	centrifuge --report-file centrifuge.report \
	-x /diskd/duh/domestic/pan_genome/centrifuge/index/nt \
	-k 1 --host-taxids 9823 \
	-f large_1000.fasta \
	> centrifuge.output
	centrifuge-kreport -x /diskd/duh/domestic/pan_genome/centrifuge/index/nt \
	centrifuge.output \
	--min-score 0 --min-length 0 \
	> centrifuge.krakenOut
done
