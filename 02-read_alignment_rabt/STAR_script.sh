#!/bin/bash


STAR limitGenomeGenerateRAM 50000000000 --runThreadN 5 --runMode genomeGenerate \
    --genomeDir ref_STAR --genomeFastaFiles Si_gnG.fa --sjdbGTFfile GCF_000188075.1_Si_gnG_genomic.gff \
    --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon exon --sjdbOverhang 115

#COMMANDS EXPLAINED
#
##STAR --runThreadN 5: Cores used for indexing
#--runMode genomeGenerate: Genome indexing mode
#--genomeDir ../ref_genome_Si_gnG/: Directory where the .fasta file for the reference genome is
##--genomeFastaFiles  path/to/reference/genome.fa : .fasta file containing reference genome
##--sjdbGTFfile path/to/annotation/file.gff3: Path to annotation file
##--sjdbGTFtagExonParentTranscript Parent: Option to be enabled if using .gff files instead of .gtf for the annotations.
#--sjdbGTFfeatureExon exon: Feature in which to group the reads. MUST BE CHANGED IF annotation file only contains CDS
##--sjdbOverhang 89: Window for generating alignements, should be = longest read length-1. In this case, average of roughly 115 -1

#Generate a text file with the name of all samples:

ls trimmed_reads/*  | cut -d "/" -f2 | sed -r  's/([A-Z]+[0-9]+_[A-Z]+_L00[0-9])(_.+)/\1/' > samples.txt

#STAR alignment, first run, for all paired samples:

parallel -j 20 'STAR --runThreadN 1 --genomeDir ref_STAR \
  --readFilesIn <(zcat split_reads/paired_{}_R1.fastq.gz) <(zcat split_reads/paired_{}_R2.fastq.gz) \
  --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFileNamePrefix bam_files/paired_{}' :::: samples.txt

#STAR alignement, second run, for all paired  samples

parallel -j 20 'STAR --runThreadN 1 --genomeDir ref_STAR \
  --readFilesIn <(zcat split_reads/paired_{}_R1.fastq.gz) <(zcat split_reads/paired_{}_R2.fastq.gz) \
  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam_files/paired_{} --outSAMstrandField intronMotif \
  --sjdbFileChrStartEnd bam_files/paired_{}SJ.out.tab' :::: samples.txt

#STAR alignment, first run, for R1 single samples:

parallel -j 20 'STAR --runThreadN 1 --genomeDir ref_STAR \
  --readFilesIn <(zcat split_reads/single_{}_R1.fastq.gz) --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
  --outFileNamePrefix bam_files/single_{}' :::: samples.txt

#STAR alignment, first run, for R2 single samples:

parallel -j 20 'STAR --runThreadN 1 --genomeDir ref_STAR \
  --readFilesIn <(zcat split_reads/single_{}_R2.fastq.gz) --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif --outFileNamePrefix bam_files/single_R2_{}' :::: samples.txt

#STAR alignement, second run, for all single R1 samples

parallel -j 20 'STAR --runThreadN 1 --genomeDir ref_STAR \
  --readFilesIn <(zcat split_reads/single_{}_R1.fastq.gz) --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix bam_files/single_{} --outSAMstrandField intronMotif \
  --sjdbFileChrStartEnd bam_files/single_{}SJ.out.tab' :::: samples.txt

#STAR alignement, second run, for all single R2 samples

parallel -j 20 'STAR --runThreadN 1 --genomeDir ref_STAR \
--readFilesIn <(zcat split_reads/single_{}_R2.fastq.gz) --outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix bam_files/single_R2_{} --outSAMstrandField intronMotif \
--sjdbFileChrStartEnd bam_files/single_{}SJ.out.tab' :::: samples.txt
