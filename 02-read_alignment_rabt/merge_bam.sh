#Merge all bam files per sample:
parallel -j 20 'samtools merge bam_files/merged_{}.bam \
  bam_files/paired_{}Aligned.sortedByCoord.out.bam bam_files/single_{}Aligned.sortedByCoord.out.bam \
  bam_files/single_R2_{}Aligned.sortedByCoord.out.bam' :::: samples.txt

#Remove all intermediate bam files generated:
rm bam_files/paired*bam
rm bam_files/single*bam
