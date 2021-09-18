#!/bin/bash
###
# @Description:
# @Author: Zhuo Yue
# @Date: 2021-07-15 23:00:43
# @LastEditors: Zhuo Yue
# @LastEditTime: 2021-08-30 11:17:43
# @Calls:
# @Called By:
# @FilePath: \CDpan\dh\5repearmasker.sh
###
#SBATCH -J repeat
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 6

for i in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S0352mu AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu; do
	cd /storage3/duh/pan_genome/Anqing/$i
	/home/duh/software/RepeatMasker/RepeatMasker -nolow -species pig /storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.final.fa
	samtools faidx /storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.final.fa
done
