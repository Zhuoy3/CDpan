#!/bin/bash
###
 # @Description:
 # @Author: Zhuo Yue
 # @Date: 2021-07-15 23:00:43
 # @LastEditors: Zhuo Yue
 # @LastEditTime: 2021-08-24 11:18:29
 # @Calls:
 # @Called By:
 # @FilePath: \CDpan\dh\4to5_1mmseqs_new.sh
###
#SBATCH -J mmseqs
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 10

source /home/duh/miniconda3/bin/activate mmseqs2
for i in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S0352mu AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu
do
	cd /storage3/duh/pan_genome/Anqing/$i/centrifuge
	time mmseqs easy-linclust \
	/storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.fa \
	/storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs \
	/storage3/duh/tmp --threads 10 --cov-mode 1 -c 0.9 --min-seq-id 0.9 --cluster-mode 2
done

centrifuge-build -p 16 --bmax 1342177280 --conversion-table gi_taxid_nucl.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp nt.fa nt
