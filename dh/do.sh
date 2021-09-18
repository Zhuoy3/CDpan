#!/bin/bash
#SBATCH -J merge
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 20

cd /storage3/duh/pan_genome/all/direct_merge/individual/534

perl /storage3/duh/pan_genome/all/individual/individual/change.pl \
/storage3/duh/pan_genome/all/direct_merge/individual/534/all.fasta \
> /storage3/duh/pan_genome/all/direct_merge/individual/534/all1.fasta

source /home/duh/miniconda3/bin/activate mmseqs2

mmseqs easy-linclust /storage3/duh/pan_genome/all/direct_merge/individual/534/all1.fasta \
/storage3/duh/pan_genome/all/direct_merge/individual/534/all1.mmseqs /storage3/duh/tmp/ \
--threads 20 --cov-mode 1 -c 0.9 --min-seq-id 0.9 --cluster-mode 2

