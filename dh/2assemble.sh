#!/bin/bash
#PBS -N test
#PSB -q coms_high
#PBS -l nodes=1:ppn=24

for i in JXH-5026 JXH-6002 JXH-7054 JXH-8062 JXH-8086 JXH-9008 JXH-9009 JXH-9028 JXH-9042 JXH-9047 JXH-9050 JXH-9057 JXH-9060 JXH-9062 JXH-9063 JXH-9068 JXH-9070 JXH-9076 JXH-9087 JXH-9088
do
	cd /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i
	touch /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/masurca_config.txt
	printf 'DATA\nPE= pe 300 50 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/'$i'/mateUnmapped_R1.fq /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/'$i'/mateUnmapped_R2.fq\nPE= s1 300 50 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/'$i'/R1_mateMapped.fq\nPE= s2 300 50 /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/'$i'/R2_mateMapped.fq\nEND\n\nPARAMETERS\nGRAPH_KMER_SIZE = auto\nUSE_LINKING_MATES = 1\nKMER_COUNT_THRESHOLD = 1\nNUM_THREADS = 24\nJF_SIZE=200000000\nDO_HOMOPOLYMER_TRIM=0\nEND' > /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/masurca_config.txt
	mkdir /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/assemble
	cd /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/assemble
	/home/liujf/software/MaSuRCA-3.3.2/bin/masurca /home/liujf/WORKSPACE/domestic_pig/pan_genome/Jiaxinghei/$i/masurca_config.txt
	./assemble.sh
done
