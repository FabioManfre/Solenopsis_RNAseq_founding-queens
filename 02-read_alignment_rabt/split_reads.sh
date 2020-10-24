#!/bin/sh

#Script for splitting into paired and single-ended from the trimmed Manfredini data:

#Unzip the trimmed reads:

parallel 'gunzip {}' ::: trimmed_reads/

#Generate a list with the names of the samples:

ls trimmed_reads/* | sed -r 's/(trimmed_reads\/)([[:alnum:]]+_[[:alnum:]]+_[[:alnum:]]+)(_.+)/\2/p' | uniq > samples.txt

#Get the IDs for paired-end reads only----------------------------------------

parallel "cat trimmed_reads/{}_R1.trimmed.fastq trimmed_reads/{}_R2.trimmed.fastq | \
  grep '^@HWI' |cut -f1 -d' '| sort |uniq -c | grep  '2 @' | cut -f2 -d '@' > split_reads/IDs_split/{}_IDs_paired.txt" :::: samples.txt

#Get the IDs for single-end reads only------------------------------------------------

parallel "cat trimmed_reads/{}_R1.trimmed.fastq trimmed_reads/{}_R2.trimmed.fastq | \
  grep '^@HWI' |cut -f1 -d' '| sort |uniq -c | grep  '1 @' | cut -f2 -d '@' > split_reads/IDs_split/{}_IDs_single.txt" :::: samples.txt


#Use seqtk for splitting the files--------------------------------------------------------

#For paired-end R1 files:

parallel 'seqtk subseq trimmed_reads/{}_R1.trimmed.fastq split_reads/IDs_split/{}_IDs_paired.txt > split_reads/paired_{}_R1.fastq' :::: samples.txt

#For paired-end R2 files:

parallel 'seqtk subseq trimmed_reads/{}_R2.trimmed.fastq split_reads/IDs_split/{}_IDs_paired.txt > split_reads/paired_{}_R2.fastq' :::: samples.txt

#gzip paired-end reads:

parallel 'gzip {}' ::: split_reads/paired_*


#For single-end R1 files:

parallel 'seqtk subseq trimmed_reads/{}_R1.trimmed.fastq split_reads/IDs_split/{}_IDs_single.txt > split_reads/single_{}_R1.fastq' :::: samples.txt

#For single-end R2 files:

parallel 'seqtk subseq trimmed_reads/{}_R2.trimmed.fastq split_reads/IDs_split/{}_IDs_single.txt > split_reads/single_{}_R2.fastq' :::: samples.txt

#Gzip single-end reads:
parallel 'gzip {}' ::: split_reads/single_*

#Gzip trimmed reads again
parallel 'gzip {}' ::: trimmed_reads/*
