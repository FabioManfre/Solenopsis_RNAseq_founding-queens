#!/bin/sh

#This generates the 'samples_names.txt' file, with only the names of the raw reads to fed them to parallel.
#This is done to allow parallel to add .trimmed.fastq directly to the trimmed reads
ls raw_reads/*gz | cut -d "/" -f2 | cut -d "." -f 1  > sample_names.txt

#This runs seqtk in 20 batches (-j 20) to remove the first 10 pb of each read. This was done accordingly to the fastqc files

parallel 'seqtk trimfq -b 10 raw_reads/{}.fastq.gz > ~/hive/SBCS-WurmLab/cmartinezruiz/sinvicta_manfredini_fastq/trimmed_reads/{}.trimmed.fastq' :::: sample_names.txt

#Remove the bad quality tiles observed in the fastqc files from the raw reads.
#This command remove the line with the matching pattern (which contains the
#information for the tiles to trim) and 3 additional lines, thus removing the whole read.

# The tiles to remove are: 2316, 2315, 2314, 2313, 2216, 2215, 2214, 2213, 2116, 2115, 2114, 2113, 1316,
#1315, 1314, 1313, 1216, 1215, 1214, 1213, 1116, 1115, 1114 and 1113. No coincidence that the worse
#looking tiles end in 16, apparently these are the tales in the borders of the flow-cell, which are more
#likely to be affected.

#fastq files need to be unzipped first

parallel 'gunzip {}' ::: trimmed_reads/*

#Sed command removing the tiles. The trimmed reads are stored in the .tmp files:
parallel -j 20 "sed -e sed -r '/@HWI-[[:alnum:]]+:[0-9]+:[[:alnum:]]+:[0-9]:(1113|1114|1115|1116|1213|1214|1215|1216|1313|1314|1315|1316|2113|2114|2115|2116|2213|2214|2215|2216|2313|2314|2315|2316):.+/,+3d' {} > {}.tmp" ::: trimmed_reads/*

#Remove the old reads:
rm trimmed_reads/*.fastq

#Rename the .tmp reads:
rename .fastq.tmp .fastq trimmed_reads/*

#Remove reads with more than 75% bp with phred score lower than 20

parallel 'fastq_quality_filter -q 20 -p 75 -z -i {} -o {}.gz' ::: trimmed_reads/*

#Remove the old reads

rm trimmed_reads/*.fastq

#And finally, remove the adapters:
#The FASTQC files showed that adapter contamination was present only at the end of the reads,
#so only 3' adapters need to be removed (-a option in cutadapt)
#The adapters were specified by the sequencing company

#Cutadapt will be run twice
#First run:
ls trimmed_reads/ | parallel -j 20 'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNATCTCGTATGCCGTCTTCTGCTTG -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed_reads/1.{} trimmed_reads/{}' :::: - &
#Remove the old reads

rm trimmed_reads/[A-Z]*[0-9]*fastq.gz

#Second run:

ls trimmed_reads/ | parallel -j 20 'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNATCTCGTATGCCGTCTTCTGCTTG -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed_reads/2.{} trimmed_reads/{}' :::: - &

#Remove reads from the first read:

rm trimmed_reads/1.*

#Rename all files

me 2\.1\. '' trimmed_reads/*
