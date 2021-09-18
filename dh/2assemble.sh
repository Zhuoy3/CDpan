#!/bin/bash
#PBS -N test
#PSB -q coms_high
#PBS -l nodes=1:ppn=24

for i in AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu AQF-S0352mu
do
	cd /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i
	touch /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/masurca_config.txt
	printf 'DATA\nPE= pe 300 50 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/'$i'/mateUnmapped_R1.fq /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/'$i'/mateUnmapped_R2.fq\nPE= s1 300 50 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/'$i'/R1_mateMapped.fq\nPE= s2 300 50 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/'$i'/R2_mateMapped.fq\nEND\n\nPARAMETERS\nGRAPH_KMER_SIZE = auto\nUSE_LINKING_MATES = 1\nKMER_COUNT_THRESHOLD = 1\nNUM_THREADS = 24\nJF_SIZE=200000000\nDO_HOMOPOLYMER_TRIM=0\nEND' > /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/masurca_config.txt
	mkdir /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/assemble
	cd /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/assemble
	/home/liujf/software/MaSuRCA-3.3.2/bin/masurca /home/liujf/WORKSPACE/domestic_pig/pan_genome/Anqing/$i/masurca_config.txt
	./assemble.sh
done
