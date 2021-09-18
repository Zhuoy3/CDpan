#!/bin/bash
#PBS -N Bamei6
#PBS -q coms_low
#PBS -l nodes=node-07:ppn=6

rawdatafile=/storage2/liujf/20190812/Bamei	##the file which rawdata locate
id=Bamei6					##the sample ID final output
rawid=Bamei6				##the ID of the raw data
line1=86			##the first line
line2=94			##the second line
output=/home/liujf/WORKSPACE/domestic_pig/Bamei/Bamei6			##the final output file

##switch to the work file
cd $output
mkdir $output/quality

##raw reads quality control
/home/liujf/software/TrimGalore-0.6.1/trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
--path_to_cutadapt /apps/.local/software/program/.virtualenvs/python_3.7.3/bin/cutadapt \
--paired $rawdatafile/$rawid\-$line1\_1.fq.gz $rawdatafile/$rawid\-$line1\_2.fq.gz \
--gzip -o $output -j 6
/home/liujf/software/TrimGalore-0.6.1/trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
--path_to_cutadapt /apps/.local/software/program/.virtualenvs/python_3.7.3/bin/cutadapt \
--paired $rawdatafile/$rawid\-$line2\_1.fq.gz $rawdatafile/$rawid\-$line2\_2.fq.gz \
--gzip -o $output -j 6
##质控参数其中-q,--length,-e以及-j可以让用户进行设置 -j为线程

##merge data
cat $output/$rawid\-$line1\_1\_val\_1.fq.gz $output/$rawid\-$line2\_1\_val\_1.fq.gz > $output/$id\_clean\_1.fq.gz
cat $output/$rawid\-$line1\_2\_val\_2.fq.gz $output/$rawid\-$line2\_2\_val\_2.fq.gz > $output/$id\_clean\_2.fq.gz

##aligning to the reference genome
/apps/bioinformatics/.local/bin/bwa mem -t 6 -M -R "@RG\tID:$id\tLB:$id\tPL:ILLUMINA\tSM:$id" \
/home/liujf/WORKSPACE/duh/11_1/reference/pig.fa \
$output/$id\_clean\_1.fq.gz $output/$id\_clean\_2.fq.gz > $output/$id.sam 
##-t为线程

##reorder sam
/apps/bioinformatics/gatk-4.0.12.0/gatk --java-options "-Xmx32G" ReorderSam \
-I $output/$id.sam -O $output/$id.reorder.sam  -R /home/liujf/WORKSPACE/duh/11_1/reference/pig.fa

##sam to bam
/apps/bioinformatics/.local/bin/samtools view -b -S \
$output/$id.reorder.sam -o $output/$id.reorder.bam

##delete process files
rm -rf $output/$rawid\-$line*
rm -rf $output/$id.reorder.sam

##sort bam
/apps/bioinformatics/gatk-4.0.12.0/gatk --java-options "-Xmx32G" SortSam \
-I $output/$id.reorder.bam -O $output/$id.sort.bam  --SORT_ORDER coordinate 
