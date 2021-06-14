# Scripts for "SOCIAL ISOLATION AND GROUP SIZE ARE ASSOCIATED WITH DIVERGENT GENE EXPRESSION IN THE BRAIN OF ANT QUEENS"

Scripts for the generation and analysis of the data used for "Complexity of the social environment and behavioural plasticity drive divergent gene expression in the brain of ant queens". The scripts listed here go through all the steps from raw RNAseq data, to read count, analyses of differences in gene expression and gene network analyses in different queen founding groups of the red fire ant *Solenopsis invicta*.

The raw data used for this analysis is available in NCBI: 
[PRJNA525584](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA525584)

All the subdirectories are numbered according to the order in which they were run for the analysis:
1. **read_qc_trimming:** Clean raw reads by quality and remove adapters.
2. **read_alignment_rabt:** Align the reads to the gnG assembly reference and generate a new RABT annotation to detect new transcripts.
3. **associate_GO_terms:** Associate GO terms from different sources to the *S. invicta* annotation.
4. **annotate_new_transcripts:** Insert the new transcripts detected in step 02 into a new annotation file.
5. **kallisto_counts:** Obtain read counts to measure expression.
6. **gene_expression_analysis:** Use DESeq2 to compare expression patterns between groups of founding queens.
7. **gene_network_analysis:** Use WGCNA to call co-expression modules in groups of founding queens.
8. **module_supergene_enrichment:** Perform an enrichment analysis to test whether any module is enriched in genes within the supergene.
9. **module_association_plots:** Perform plots to check how different genes change across modules.
