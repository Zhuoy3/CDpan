#!/bin/bash

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088 
do
	if test -s /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.delta.coords; then
		tail -n +6 /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.delta.coords | \
		awk '{print $19}' - | sort -u - > /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/repeat.names
		fasta_tool --remove /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/repeat.names \
		/storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs_rep_seq.fasta \
		> /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.final.fa
	else
		cp /storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs_rep_seq.fasta \
		/storage3/duh/pan_genome/Jiaxinghei/$i/centrifuge/filtered.mmseqs.final.fa
	fi
done
