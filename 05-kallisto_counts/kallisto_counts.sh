#!/bin/sh

#Generate a Kallisto index for the SignG annotation as it is from NCBI:

kallisto index -i k_index ref_gnG/s_invicta_gnG_transcriptome.fa

#Generate a list with the names of the samples:
ls trimmed_reads/*fastq.gz | sed -r 's/(trimmed_reads\/)([[:alnum:]]+_[[:alnum:]]+_[[:alnum:]]+)(_.+)/\2/p'|uniq > samples.txt

parallel 'kallisto quant -i k_index_manfredini -o ./{} --plaintext trimmed_reads/{}_R1.trimmed.fastq.gz trimmed_reads/{}_R2.trimmed.fastq.gz' :::: samples.txt

#-------------------------- Using the RABT annotation file from cufflinks --------------------------------------

#Generate a Kallisto index, generated in step 02:

kallisto index -i k_index_cuff ref_gnG/s_invicta_gnG_transcriptome_cuff.fa

#Generate a list with the names of the samples:
ls trimmed_reads/*fastq.gz | sed -r 's/(trimmed_reads\/)([[:alnum:]]+_[[:alnum:]]+_[[:alnum:]]+)(_.+)/\2/p'|uniq > samples.txt

parallel 'kallisto quant -i k_index_manfredini_cuff -o ./{}_cuff --plaintext trimmed_reads/{}_R1.trimmed.fastq.gz trimmed_reads/{}_R2.trimmed.fastq.gz' :::: samples.txt
