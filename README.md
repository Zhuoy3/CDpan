# CDpan: Constructing dispensable genome for pan-genome analysis

## Introduction

There is an increasing understanding that a reference sequence representing a genome of individual is insufficient to capture the genomic diversity we observe in nature. With the development of short-reads sequencing, it’s convenient to achieve the high depth whole-genome sequencing data of an individual. Therefore, to help the researchers to capture the genomic sequences absent in the reference genome within the species, we propose the CDpan (Constructing dispensable genome for pan-genome analysis) to construct the dispensable genome using high-throughput sequencing reads.

## Installation

Download from project homepage, <https://github.com/kimi-du-bio/CDpan>, or check out the code using Git from <https://github.com/kimi-du-bio/CDpan.git>.

### Dependencies

CDpan require Perl version 5.16 or greater to run. Only Linux is supported (May or may not run under Perl for MacOS or Windows, etc).

CDpan uses three Perl modules, ```Bio::SeqIO```, ```Config::IniFiles``` and ```File::Slurp```, which may be not installed by default configuration of perl, you can install them as following:

    cpan -i Bio::SeqIO Config::IniFiles File::Slurp

Using CDpan requires the following tools to be installed:

- [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
- [Cutadapt](https://github.com/marcelm/cutadapt)
- [FastQC](https://github.com/s-andrews/FastQC)
- [BWA]()
- [GATK]()
- [samtools]()
- [MaSuRCA 3.x]() (Not compatible with MaSuRCA 4.0 or greater)
- [centrifuge]()
- [MMseqs2]()
- [mummer]()
- [RepeatMasker]()
- [Bowtie2]()
- [bedtools]()
- [minimap2]()

CDpan has been tested on the following environment: ```Perl v5.34, TrimGalore v0.6.7, Cutadapt v3.4, FastQC v0.11.9, BWA v0.7.17, GATK v4.2.1.0, samtools v1.13, MaSuRCA v3.4.2, centrifuge v1.0.4, MMseqs2 v13, mummer v4.0.0, RepeatMasker v4.1.2, Bowtie2 v2.4.4, bedtools v2.30.0, minimap2 v2.22```. It may run on any other version of most tools except the MaSuRCA 4.0 or greater.

Building the latest version from the repository
Quick Start Guide
Running CDpan
Run CDpan from the command line like this:

    CDpan <module> -i input_path -c config_file [options]
    Module: filter
            align
            extract
            assembly
            mope
            vot
            soot
            merge
            location

            RUN-ALL            Run all modules of CDpan
            RUN-DIEM       Run all modules of CDpan except location

    Options: -i, --input             path of input file (Required)
             -c, --config            path of config file (Required)
             -w, --work_dir          path of work directory ( Default is directory \'./cdpan_tmp\' )
             -p, --process           whether to keep process files (Default is False)
             -o, --output                prefix of the output file (Default is prefix of input file)
             -O, --output_dir        output directory (Default is the current directory)
                 --no-qc             no quality control (Default is False)
             -v, --version           print version message
             -h, --help              print help message

The main input of the program is the fastq file with paired-end sequencing reads.
The main output of the program is a fasta file containing the sequences of the constructed dispensable genome and the predicted inserted location of the new sequences relative to the reference genome.
To use CDpan, you need to specify a module and support the "input_path". Please refer to the "Usage" section to get the detailed information of the module contained in the CDpan. The introduction of "input_path" please refer to the following "Usage" section.

3.Usage
3.1 CDpan quick usage

    CDpan RUN-ALL -i input_file -c config_file

Input_path
The "input_path" is the directory which you archive your sequencing data. The sequencing data must be the paired-end reads, the directory format should be like this:
Data
|
|---sample1---|---*_1.fastq (fq, fastq.gz, fq.gz)
              |---*_2.fastq (fq, fastq.gz, fq.gz)
|---sample2---|---*_1.fastq (fq, fastq.gz, fq.gz)
              |---*_2.fastq (fq, fastq.gz, fq.gz)
……
|---sampleN---|---*_1.fastq (fq, fastq.gz, fq.gz)
              |---*_2.fastq (fq, fastq.gz, fq.gz)

config_file
The "config_file" can be found in the XXX folder, and change this file to suit your practical circumstance. The detailed described of the "config_file", please refer to the "Config file" section.

Output
This program produces a fasta file containing the sequences of the constructed dispensable genome and the txt file for the predicted inserted location.

If you only interested in the DNA sequences absent in the reference genome but present within the species, you can use the following command:

    CDpan RUN-DIEM -i input_file -c config_file

The output is only the fasta file containing the new DNA sequences.

3.2 CDpan (traditional) usage
Flexible and customizable new DNA sequences detection for advanced users.

Usage:
    CDpan <module> -i input_path -c config_file [options]
    Module: filter
            align
            extract
            assembly
            mope
            vot
            soot
            merge
            location

    Options: -i, --input         path of input file (Required)
             -c, --config        path of config file (Required)
             -w, --work_dir      path of work directory ( Default is directory \'./cdpan_tmp\' )
             -p, --process       whether to keep process files (Default is False)
             -o, --output        prefix of the output file (Default is prefix of input file)
             -O, --output_dir    output directory (Default is the current directory)
                 --no-qc         no quality control (Default is False)
             -v, --version       print version message
             -h, --help          print help message

Note:Before running the follwong steps, we recommend to first set the "config_file". The following steps use the same "config_file".The detailed information of the "config_file" can be found in the "Config file" section.

##Pre-processing

    CDpan filter -i input_path -c config_file

The "input_path" is the same format as the quick usage. This module is used to filter the low quality reads. You can do the quality control using other software and run the following processes.

##Alignment

    CDpan align -i input_path -c config_file

The "input_path" is the same format as the quick usage. If you use the CDpan to do the pre-processing procedure, the "input_path" is the output directory of the pre-processing. Of course, you can use the clean data which processed by yourself, just provide the data deposited directory to the "input_path". This module is used to align the sequencing reads to the reference genome.

##Extract the unmapped reads

    CDpan extract -i input_path -c config_file

The "input_path" is the output directory of the alignment process. This step is used to extract the ummaped reads from the aligned results.

##Assembly

    CDpan assembly -i input_path -c config_file

The "input_path" is the output directory of the previous process. This step is used to assemble the unmmaped reads of each individual.

##Remove the contaminants

    CDpan mope -i input_path -c config_file

The "input_path" is the output directory of the assembly process. This step is used to remove the contaminants sequences in the assembled results.

##Remove the redundants

    CDpan vot -i input_path -c config_file

The "input_path" is the output directory after the "mpope" module running. This step is used to remove the redundants sequences in the assembled results.

##Retained the new DNA sequences

    CDpan soot -i input_path -c config_file

The "input_path" is the output directory of the "vot" process. This step is used to keep the absent sequences of the reference genome in the assembled results.

##Final dispensable genome

    CDpan merge -i input_path -c config_file

The "input_path" is the output directory of the "soot" process. This step is used to integrated the new DNA sequences of all individuals and construct the final dispensable genome. The output is a fasta file.

##Location

    CDpan location -i input_path -c config_file

The "input_path" is XXXX. This step is used to determine the genomic positions of new DNA sequences relative to the reference genome. The output of this step is a txt file.

4.Config file
The follwing is a example of the "config_file", you can also find this file in the XXX folder:

# This part definites the threads used in the pipeline

[CDPAN]
thread = 1

##This part definites the path of the dependency software
[TOOLS]
trim_galore = null
cutadapt = null
fastqc = null
bwa = null
gatk = null
samtools = null
masurca = null
centrifuge = null
centrifuge-kreport = null
mmseqs = null
nucmer = null
show-coords = null
RepeatMasker = null
bowtie2 = null
bedtools = null
bowtie2-build = null

##This section definites the required file
[DATA]
ref = *.fasta       ##The reference genome file, must be .fasta or .fa format
nt_index =          ##The NCBI nucleotide non-redundant index, download from https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
taxid =             ##The file contained the NCBI taxid which you want to keep, just one colunm. We supply a example file contained all chordate taxid in the XXX.

##The parameter of pre-processing.
[FILTER]
quality  = 20       ##The minimum Phred score
length = 20         ##The minimum length of reads after quality and adapter trimming
error-rate = 0.1    ##Maximum allowed error rate


[ALIGN]

[EXTRACT]

##The parameter of assembly.
[ASSEMBLY]
fragment-mean=300   ##The fragment mean of paired-end reads
fragment-stdev=50   ##The fragment standard deviation of paired-end reads
JF_SIZE=jellyfish hash size, set this to about 10x the genome size.

##The parameter of mope module.
[MOPE]
host-taxids = null  ##The NCBI taxid of your research species
min-length = 1000   ##The minmum length of new DNA sequences

##The parameter of vot module (we recommend to only just change the "min-seq-id" parameter).
[VOT]
cov-mode = 1
coverage = 0.9
min-seq-id = 0.9    ##The ratio of length overlap fortwo sequences
cluster-mode = 2

##The parameter of soot module
[SOOT]
maxgap = 1000       ##The maximum gap between two adjacent matches in a cluster
mincluster = 90     ##The minimum length of a cluster of matches
minmatch = 500      ##The minimum length of a single exact match

[MERGE]

[LOCATION]
species = null

<!--
Features to be added
1.Use perl language to rewrite ex.py
2.Log
3.Almost all parameters could be imported using the parameter card
4.Safer file path handling
5.Allows to import reference genomes that already have dictionary and index
-->
