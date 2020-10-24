library(tximport)
library(readr)

########################################################################LOAD KALLISTO DATA INTO R USING Tximport#######################################################################

# Load Kallisto abundance files---------------------------------------------------------------------------
kallisto_files <- grep(".tsv", list.files("Path/to/Kallisto_counts", recursive = TRUE), value = TRUE)
kallisto_files <- paste("Path/to/Kallisto_counts/", kallisto_files, sep="" )

#Extract the sample names:
kallisto_names <- gsub(pattern = "(Path/to/Kallisto_counts/)([[:alnum:]]+)(_[[:alnum:]]+_[[:alnum:]]+)(\\\\.+)",
                       "\\2", kallisto_files)

#Label the files:
names(kallisto_files) <- kallisto_names

#Generate the colData table for DESeq2 analysis (information given by Fabio):
data_info <- read.csv("Path/to/samples_treatment.csv", sep = "\t", header = TRUE)

#Generate an additional column for merged factors (Haplometrotic + Pleometric -A2+A3- and Pleometrotic big + Pleometrotic small -B3 + B4-)

library(rockchalk)
data_info$merged <- combineLevels(data_info$treatment,levs = c("A2", "A3", "B3", "B4"), newLabel = c("grouped"))

#Separate the kallisto files into the two different experiments:

kallisto_files_A <- kallisto_files[names(kallisto_files) %in% data_info$sample_name[grep("A", data_info$treatment)]]
kallisto_files_B <- kallisto_files[names(kallisto_files) %in% data_info$sample_name[grep("B", data_info$treatment)]]

#Separate the data in the data_info

data_info_A <- data_info[grep("A", data_info$treatment), ]
data_info_B <- data_info[grep("B", data_info$treatment), ]

#Load a table with the transcript ID and the gene to which they belong.
#This file was generated in the directory /data/home/btx076/2016-05-09_DC_S.invicta_project/2016-11-29_rna_id_gnG in Apocrita
gene_ID <- read.table("Path/to/transcripts_and_genes_mit.txt", fill = TRUE)
colnames(gene_ID)<-c("Transcript", "Gene")
gene_ID$Transcript <- as.character(gene_ID$Transcript)
gene_ID$Gene <- as.character(gene_ID$Gene)

#Gene ID has a blank space for those transcripts that do not have an associated gene.
#These blanks are replaced by the name of the transcript
for (row in 1:nrow(gene_ID)){
  if (gene_ID[row, 2]==""){
    gene_ID[row, 2] <- gene_ID[row, 1]
  }
}
#The first elements from the "Transcript" column are actually mitochondrial genes, which were assigned differently in the gff file. 

#Import the abundance files into a DESeq2 object:
txi_reads_A <- tximport(kallisto_files_A, type="kallisto", tx2gene = gene_ID, reader=read_tsv)
txi_reads_B <- tximport(kallisto_files_B, type="kallisto", tx2gene = gene_ID, reader=read_tsv)

#txi_reads$abundance retrieves the TPMs for all samples


##############################################################################DESeq2 ANALYSIS####################################################################
library(DESeq2)
#Load into a DESeq2 object:
ddsTxi_A <- DESeqDataSetFromTximport(txi_reads_A, colData= data_info_A, design = ~treatment)
ddsTxi_B <- DESeqDataSetFromTximport(txi_reads_B, colData= data_info_B, design = ~treatment)
#Performs DESeq2 test for detecting DE loci, using first a likelihood ratio test (ANOVA-like). This test compares the model ~Treatment against a
#reduced model, ~1. The p values indicate which genes are significantly DE accross ALL levels of Treatment----------------------------------------------------------------------------------------- 

DE_txi_LRT_A <- DESeq(ddsTxi_A, test = "LRT", reduced = ~1)
DE_txi_LRT_B <- DESeq(ddsTxi_B, test = "LRT", reduced = ~1)
#Get the results from the dds object:
resDE_txi_LRT_A <- results(DE_txi_LRT_A)
resDE_txi_LRT_B <- results(DE_txi_LRT_B)

#The p values and corrected p values obtained here are those for the whole comparison, only the LFC values are comparison-specific. In this case, treatment B4 vs A.B1.
#In short, p values are ok, LFC values can be ignored here.

#Output the results for LRT-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Generate the dataframe for output:

LRT_p_val_A <- data.frame(resDE_txi_LRT_A$pvalue, resDE_txi_LRT_A$padj, row.names = rownames(resDE_txi_LRT_A))
colnames(LRT_p_val_A) <- c("LRT_p_values_experiment_A", "LRT_adj_p_values_experiment_A")

LRT_p_val_B <- data.frame(resDE_txi_LRT_B$pvalue, resDE_txi_LRT_B$padj, row.names = rownames(resDE_txi_LRT_B))
colnames(LRT_p_val_B) <- c("LRT_p_values_experiment_B", "LRT_adj_p_values_experiment_B")

#Outout the .txt files
write.table(LRT_p_val_A, "Path/to/output/Results_LRT_A.txt")
write.table(LRT_p_val_B, "Path/to/Results_LRT_B.txt")

#Pairwise comparisons-------------------------------------------------------------------------------------------------------------------------

DE_txi_A <- DESeq(ddsTxi_A)
DE_txi_B <- DESeq(ddsTxi_B)

#Pairwise comparisons for A:
res_DE_A_1_v_2 <- results(DE_txi_A, contrast = c("treatment", "A/B1", "A2"))
res_DE_A_1_v_3 <- results(DE_txi_A, contrast = c("treatment", "A/B1", "A3"))
res_DE_A_2_v_3 <- results(DE_txi_A, contrast = c("treatment", "A2", "A3"))

#Pairwise comparison for A2 and A3 (Haplo and Pleometrotic) combined together:
ddsTxi_A_grouped <- DESeqDataSetFromTximport(txi_reads_A, colData= data_info_A, design = ~merged)
DE_txi_A_grouped <- DESeq(ddsTxi_A_grouped)
res_DE_A_1_v_2_3 <- results(DE_txi_A_grouped)

#Pairwise comparisons for B:
res_DE_B_1_v_2 <- results(DE_txi_B, contrast = c("treatment", "A/B1", "B2"))
res_DE_B_1_v_3 <- results(DE_txi_B, contrast = c("treatment", "A/B1", "B3"))
res_DE_B_1_v_4 <- results(DE_txi_B, contrast = c("treatment", "A/B1", "B4"))
res_DE_B_3_v_4 <- results(DE_txi_B, contrast = c("treatment", "B3", "B4"))
res_DE_B_2_v_3 <- results(DE_txi_B, contrast = c("treatment", "B2", "B3"))
res_DE_B_2_v_4 <- results(DE_txi_B, contrast = c("treatment", "B2", "B4"))

#Pairwise comparisons for A3 and A4 (Pleometrotic big and small) merged together
ddsTxi_B_grouped <- DESeqDataSetFromTximport(txi_reads_B, colData= data_info_B, design = ~merged)
DE_txi_B_grouped <- DESeq(ddsTxi_B_grouped)

res_DE_B_1_v_3_4 <- results(DE_txi_B_grouped, contrast = c("merged", "A/B1", "grouped"))
res_DE_B_2_v_3_4 <- results(DE_txi_B_grouped, contrast = c("merged", "B2", "grouped"))

#Output the results for comparisons------------------------------------------------------------------------------------------------------------

#For experiment A
results_pairwise_A <- data.frame(res_DE_A_1_v_2$log2FoldChange, res_DE_A_1_v_2$pvalue, res_DE_A_1_v_2$padj, res_DE_A_1_v_3$log2FoldChange, 
                                 res_DE_A_1_v_3$pvalue, res_DE_A_1_v_3$padj, res_DE_A_2_v_3$log2FoldChange, res_DE_A_2_v_3$pvalue, res_DE_A_2_v_3$padj, 
                                 res_DE_A_1_v_2_3$log2FoldChange, res_DE_A_1_v_2_3$pvalue, res_DE_A_1_v_2_3$padj, row.names = rownames(DE_txi_A))

#Change column names:

colnames(results_pairwise_A) <- paste(rep(c("Time0_Haplo", "Time0_Pleo", "Haplo_Pleo", "Time0_Haplo+Pleo"), each = 3),
                                      rep(c("LFC", "p_value", "p_adj"), 3), sep = "_")

#For experiment B

results_pairwise_B <- data.frame(res_DE_B_1_v_2$log2FoldChange, res_DE_B_1_v_2$pvalue, res_DE_B_1_v_2$padj, res_DE_B_1_v_3$log2FoldChange, 
                                 res_DE_B_1_v_3$pvalue, res_DE_B_1_v_3$padj, res_DE_B_1_v_4$log2FoldChange, res_DE_B_1_v_4$pvalue, res_DE_B_1_v_4$padj,
                                 res_DE_B_3_v_4$log2FoldChange, res_DE_B_3_v_4$pvalue, res_DE_B_3_v_4$padj, res_DE_B_2_v_3$log2FoldChange, 
                                 res_DE_B_2_v_3$pvalue, res_DE_B_2_v_3$padj, res_DE_B_2_v_4$log2FoldChange, res_DE_B_2_v_4$pvalue, res_DE_B_2_v_4$padj,
                                 res_DE_B_1_v_3_4$log2FoldChange, res_DE_B_1_v_3_4$pvalue, res_DE_B_1_v_3_4$padj, res_DE_B_2_v_3_4$log2FoldChange, 
                                 res_DE_B_2_v_3_4$pvalue, res_DE_B_2_v_3_4$padj, row.names = rownames(DE_txi_B))

colnames(results_pairwise_B) <- paste(rep(c("Time0_Haplo", "Time0_Pleo_small", "Time0_Pleo_big", "Pleo_big_Pleo_small", "Haplo_Pleo_small", 
                                            "Haplo_Pleo_big", "Time0_Pleo", "Haplo_Pleo"), each = 3),
                                      rep(c("LFC", "p_value", "p_adj"), 3), sep = "_")

#Outpout the .txt files
write.table(results_pairwise_A, "Path/to/Results_pairwise_A.txt")
write.table(results_pairwise_B, "Path/to/Results_pairwise_B.txt")

#Output table with raw counts, normalised counts by DESeq2 and TPM per gene per sample:------------------------------------------------------------------

#Generate a vector with all the raw counts:
raw_counts_A <- as.vector(counts(DE_txi_A, normalized = FALSE))
raw_counts_B <- as.vector(counts(DE_txi_B, normalized = FALSE))

#Generate a vector with all normalised counts: 
norm_counts_A <- as.vector(counts(DE_txi_A, normalized = TRUE))
norm_counts_B <- as.vector(counts(DE_txi_B, normalized = TRUE))

#Generate a vector with all TPMs:
TPMs_A <- as.vector(txi_reads_A$abundance)
TPMs_B <- as.vector(txi_reads_B$abundance)

#Generate a vector with the sample names:
samples_A <- rep(colnames(txi_reads_A$abundance), each = nrow(txi_reads_A$abundance))
samples_B <- rep(colnames(txi_reads_B$abundance), each = nrow(txi_reads_B$abundance))

#Generate a vector with the gene names:
genes_A <- rep(rownames(txi_reads_A$abundance), ncol(txi_reads_A$abundance))
genes_B <- rep(rownames(txi_reads_B$abundance), ncol(txi_reads_B$abundance))

#Generate a vector with the levels per sample (in DESeq everything is sorted alphabetically per sample, so):
sorted_data_info_A <- data_info_A[order(data_info_A$sample_name), ]
levels_A <- rep(sorted_data_info_A$treatment, each = nrow(txi_reads_A$abundance))

sorted_data_info_B <- data_info_B[order(data_info_B$sample_name), ]
levels_B <- rep(sorted_data_info_B$treatment, each = nrow(txi_reads_B$abundance))

#Generate the actual data frame
data_experiment_A <- data.frame(genes_A, samples_A, levels_A, raw_counts_A, norm_counts_A, TPMs_A)
data_experiment_B <- data.frame(genes_B, samples_B, levels_B, raw_counts_B, norm_counts_B, TPMs_B)

#Outpout the .txt files
write.table(data_experiment_A, "Path/to/data_experiment_A.txt")
write.table(data_experiment_B, "Path/to/data_experiment_B.txt")

#OBP17.2 is significant by very little (p value = 0.049 using obp15 transcriptome), so I guess any modification can make it non-significant even if the counts are the same
#OBP17.2 will not be counted as significant