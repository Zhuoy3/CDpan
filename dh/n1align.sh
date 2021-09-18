#!/bin/bash
#SBATCH -J align
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=50gb

for id in AQ-1 AQ-2 AQ-A117mu AQ-D037mu AQF-B0079gong AQF-B0088gong AQF-B0090gong AQF-B0092gong AQF-S0352mu AQF-S01163mu AQF-S01164mu AQF-S01172 AQF-S01174mu AQF-S01180mu AQF-S01185mu
do
	mkdir /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id
	cd /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id
	bwa mem -t 10 -M -R "@RG\tID:$id\tLB:$id\tPL:ILLUMINA\tSM:$id" \
	/storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/dispensable_genome.fasta \
	/storage3/duh/pan_genome/Anqing/$id/unaligned_reads/mateUnmapped_R1.fq \
	/storage3/duh/pan_genome/Anqing/$id/unaligned_reads/mateUnmapped_R2.fq \
	> /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.sam 

	##reorder sam
	/home/duh/software/gatk-4.1.2.0/gatk --java-options "-Xmx32G" ReorderSam \
	-I /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.sam \
	-O /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.reorder.sam \
	-R /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/dispensable_genome.fasta

	##sam to bam
	samtools view -b -S \
	/storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.reorder.sam \
	-o /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.reorder.bam

	##delete process files
	rm -rf /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.reorder.sam

	##sort bam
	/home/duh/software/gatk-4.1.2.0/gatk --java-options "-Xmx32G" SortSam \
	-I /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.reorder.bam \
	-O /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.sort.bam \
	--SORT_ORDER coordinate 

	##delete process file
	rm -rf /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.reorder.sam
	rm -rf /storage3/duh/pan_genome/all/direct_merge/final/annotation/20210525/results/11_1/uniprot/align/Anqing/$id/$id.sam
done
