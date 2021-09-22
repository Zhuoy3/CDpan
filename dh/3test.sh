#!/bin/bash
#PBS -N Jiaxinghei
#PBs -q batch

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088
do
	cd /diskd/duh/domestic/pan_genome/Jiaxinghei/$i
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