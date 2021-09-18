#!/bin/bash
###
# @Description:
# @Author: Zhuo Yue
# @Date: 2021-07-15 23:00:43
# @LastEditors: Zhuo Yue
# @LastEditTime: 2021-08-31 17:45:29
# @Calls:
# @Called By:
# @FilePath: \CDpan\dh\4to5_3_de_repeat.sh
###

for i in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S0352mu AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu; do
	if test -s /storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.delta.coords; then
		tail -n +6 /storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.delta.coords |
			awk '{print $19}' - | sort -u - >/storage3/duh/pan_genome/Anqing/$i/centrifuge/repeat.names
		fasta_tool --remove /storage3/duh/pan_genome/Anqing/$i/centrifuge/repeat.names \
			/storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs_rep_seq.fasta \
			>/storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.final.fa
	else
		cp /storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs_rep_seq.fasta \
			/storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.final.fa
	fi
done
