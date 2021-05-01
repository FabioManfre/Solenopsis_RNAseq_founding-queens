library(readxl)
library(tximport)
library(readr)
library(DESeq2)
library(ggbiplot)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(psych)
library(reshape2)
library(ggbeeswarm)
library(rockchalk)

########################################################################LOAD KALLISTO DATA INTO R USING Tximport#######################################################################

# Load Kallisto abundance files---------------------------------------------------------------------------
kallisto_files <- grep(".tsv", list.files("Path/to/Kallisto_counts", recursive = TRUE), value = TRUE)
kallisto_files <- paste("Path/to/Kallisto_counts/", kallisto_files, sep="" )

#Extract the sample names:
kallisto_names <- gsub(pattern = "(Path/to/Kallisto_counts/)([[:alnum:]]+)(_[[:alnum:]]+_[[:alnum:]]+)(\\\\.+)",
                       "\\2", kallisto_files)

#Label the files:
names(kallisto_files) <- kallisto_names

#Generate the colData table for DESeq2:
all_metadata <- read_xlsx("path/to/samples_annotations_master-file.xlsx")

data_info <- all_metadata[c("Sample Name", "Treatment", "TRAY ID")]
data_info <- data.frame(data_info)
colnames(data_info) <- c("sample_name", "treatment", "tray_id")
rownames(data_info) <- data_info$sample_name

#Rename treatment according to publication
data_info$treatment_raw <- data_info$treatment
data_info$treatment <- mapvalues(data_info$treatment_raw,
                                 from = c("A/B1", "A2", "A3", "B2", "B3", "B4"),
                                 to = c("NMQ", "SFQ_3dpf", "GFQ_3dpf", "SFQ_25dpf", "GFQsmall_25dpf", "GFQlarge_25dpf"))

#Generate an additional column for merged factors (Haplometrotic + Pleometric -A2+A3- and Pleometrotic big + Pleometrotic small -B3 + B4-)
data_info$merged <- combineLevels(data_info$treatment_raw, levs = c("A2", "A3", "B3", "B4"), newLabel = c("grouped"))

#Separate the kallisto files into the two different experiments:

kallisto_files_A <- kallisto_files[names(kallisto_files) %in% data_info$sample_name[grep("A", data_info$treatment_raw)]]
kallisto_files_B <- kallisto_files[names(kallisto_files) %in% data_info$sample_name[grep("B", data_info$treatment_raw)]]

#Separate the data in the data_info

data_info_A <- data_info[grep("A", data_info$treatment_raw), ]
data_info_A <- data_info_A[order(data_info_A$sample_name), ]
data_info_B <- data_info[grep("B", data_info$treatment_raw), ]
data_info_B <- data_info_B[order(data_info_B$sample_name), ]

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

##############################################################################DESeq2 ANALYSIS####################################################################
#Load into a DESeq2 object:
ddsTxi_A <- DESeqDataSetFromTximport(txi_reads_A, colData= data_info_A, design = ~treatment_raw)
ddsTxi_B <- DESeqDataSetFromTximport(txi_reads_B, colData= data_info_B, design = ~treatment_raw)
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

#Pairwise comparisons-------------------------------------------------------------------------------------------------------------------------

DE_txi_A <- DESeq(ddsTxi_A)
DE_txi_B <- DESeq(ddsTxi_B)

#Pairwise comparisons for A:
res_DE_A_1_v_2 <- results(DE_txi_A, contrast = c("treatment_raw", "A/B1", "A2"))
res_DE_A_1_v_3 <- results(DE_txi_A, contrast = c("treatment_raw", "A/B1", "A3"))
res_DE_A_2_v_3 <- results(DE_txi_A, contrast = c("treatment_raw", "A2", "A3"))

#Pairwise comparison for A2 and A3 (Haplo and Pleometrotic) combined together:
ddsTxi_A_grouped <- DESeqDataSetFromTximport(txi_reads_A, colData= data_info_A, design = ~merged)
DE_txi_A_grouped <- DESeq(ddsTxi_A_grouped)
res_DE_A_1_v_2_3 <- results(DE_txi_A_grouped)

#Pairwise comparisons for B:
res_DE_B_1_v_2 <- results(DE_txi_B, contrast = c("treatment_raw", "A/B1", "B2"))
res_DE_B_1_v_3 <- results(DE_txi_B, contrast = c("treatment_raw", "A/B1", "B3"))
res_DE_B_1_v_4 <- results(DE_txi_B, contrast = c("treatment_raw", "A/B1", "B4"))
res_DE_B_3_v_4 <- results(DE_txi_B, contrast = c("treatment_raw", "B3", "B4"))
res_DE_B_2_v_3 <- results(DE_txi_B, contrast = c("treatment_raw", "B2", "B3"))
res_DE_B_2_v_4 <- results(DE_txi_B, contrast = c("treatment_raw", "B2", "B4"))

#Pairwise comparisons for A3 and A4 (Pleometrotic big and small) merged together
ddsTxi_B_grouped <- DESeqDataSetFromTximport(txi_reads_B, colData= data_info_B, design = ~merged)
DE_txi_B_grouped <- DESeq(ddsTxi_B_grouped)

res_DE_B_1_v_3_4 <- results(DE_txi_B_grouped, contrast = c("merged", "NMQ", "grouped"))
res_DE_B_2_v_3_4 <- results(DE_txi_B_grouped, contrast = c("merged", "SFQ_25dpf", "grouped"))

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

#Perform PCA on the expression data---------------------------------------------------------------------------------------------------
vsd_A <- vst(DE_txi_A, blind = TRUE)
vsd_B <- vst(DE_txi_B, blind = TRUE)

#Extract VST counts
vst_counts_A <- assay(vsd_A)
vst_counts_B <- assay(vsd_B)

#Get KMO values
kmo_A <- KMO(vst_counts_A)
kmo_B <- KMO(vst_counts_B)

#Perform PCA
pdf(file = "results/pca_varA.pdf", width = 10, height = 7.5)
pca_A <- prcomp(t(vst_counts_A))
plot_pca_A <- stats:::summary.prcomp(pca_A)$importance[2, ]
barplot(plot_pca_A, ylab = "Proportion of variance",
        ylim = c(0, 0.5), xlim = c(0, 21.7), las = 2)
dev.off()

#Save all results as tables
pcaA_loadings <- pca_A$rotation
pcaA_eigenvalues <- pca_A$x

write.csv(pcaA_loadings, file = "results/pcaA_loadings.csv", quote = FALSE)
write.csv(pcaA_eigenvalues, file = "results/pcaA_eigenvalues.csv", quote = FALSE)

pdf(file = "results/pca_varB.pdf", width = 10, height = 7.5)
pca_B <- prcomp(t(vst_counts_B))
plot_pca_B <- stats:::summary.prcomp(pca_B)$importance[2, ]
barplot(plot_pca_B, ylab = "Proportion of variance",
        ylim = c(0, 0.5), xlim = c(0, 21.7), las = 2)
dev.off()

#Save all results as tables
pcaB_loadings <- pca_B$rotation
pcaB_eigenvalues <- pca_B$x

write.csv(pcaB_loadings, file = "results/pcaB_loadings.csv", quote = FALSE)
write.csv(pcaB_eigenvalues, file = "results/pcaB_eigenvalues.csv", quote = FALSE)

#Plot by treatment and tray
pdf(file = "results/PCA_A_treatment.pdf", width = 10, height = 8)
ggbiplot(pca_A, obs.scale = 1, var.scale = 1,
         groups = data_info_A$treatment, ellipse = TRUE, var.axes = FALSE,
         circle = TRUE) + scale_shape_discrete(name = '') +
  geom_point((aes(color = data_info_A$treatment)), size = 2) +
  theme_bw() + scale_color_manual(values = brewer.pal(3, "Set1"),
                                  name = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

pdf(file = "results/PCA_A_tray.pdf", width = 10, height = 8)
ggbiplot(pca_A, obs.scale = 1, var.scale = 1,
         groups = data_info_A$tray_id, ellipse = FALSE, var.axes = FALSE,
         circle = FALSE) + scale_shape_discrete(name = '') +
  geom_point((aes(color = data_info_A$tray_id)), size = 2) +
  theme_bw() + scale_color_discrete(name = "Tray ID") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

pdf(file = "results/PCA_B_treatment.pdf", width = 10, height = 8)
ggbiplot(pca_B, obs.scale = 1, var.scale = 1,
         groups = data_info_B$treatment, ellipse = TRUE, var.axes = FALSE,
         circle = TRUE) + scale_shape_discrete(name = '') +
  geom_point((aes(color = data_info_B$treatment)), size = 2) +
  theme_bw() + scale_color_manual(values = brewer.pal(4, "Set2"),
                                  name = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

pdf(file = "results/PCA_B_tray.pdf", width = 10, height = 8)
ggbiplot(pca_B, obs.scale = 1, var.scale = 1,
         groups = data_info_B$tray_id, ellipse = FALSE, var.axes = FALSE,
         circle = FALSE) + scale_shape_discrete(name = '') +
  geom_point((aes(color = data_info_B$tray_id)), size = 2) +
  theme_bw() + scale_color_discrete(name = "Tray ID") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

#Repeat PCA removing NA tray samples (control)--------------------------------------------------------------------------------------------
na_tray_samples <- data_info$sample_name[data_info$tray_id == "n.a."]
data_info_A_nona <- data_info_A[!data_info_A$sample_name %in% na_tray_samples, ]
vst_counts_A_nona_tray <- vst_counts_A[, !colnames(vst_counts_A) %in% na_tray_samples]

pca_A_nona <- prcomp(t(vst_counts_A_nona_tray))
plot_pca_A_nona <- stats:::summary.prcomp(pca_A_nona)$importance[2, ]
barplot(plot_pca_A_nona, ylab = "Proportion of variance",
        ylim = c(0, 0.25), xlim = c(0, 21.7), las = 2)


pdf(file = "results/PCA_A_treat_nona.pdf", width = 10, height = 8)
ggbiplot(pca_A_nona, obs.scale = 1, var.scale = 1,
         groups = data_info_A_nona$treatment, ellipse = FALSE, var.axes = FALSE,
         circle = FALSE) + scale_shape_discrete(name = '') +
  geom_point((aes(color = data_info_A_nona$treatment)), size = 2) +
  theme_bw() + scale_color_manual(values = brewer.pal(3, "Set1"),
                                  name = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

pdf(file = "results/PCA_A_tray_nona.pdf", width = 10, height = 8)
ggbiplot(pca_A_nona, obs.scale = 1, var.scale = 1,
         groups = data_info_A_nona$tray_id, ellipse = FALSE, var.axes = FALSE,
         circle = FALSE) + scale_shape_discrete(name = '') +
  geom_point((aes(color = data_info_A_nona$tray_id)), size = 2) +
  theme_bw() + scale_color_discrete(name = "Tray ID") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

#Same with B
data_info_B_nona <- data_info_B[!data_info_B$sample_name %in% na_tray_samples, ]
vst_counts_B_nona_tray <- vst_counts_B[, !colnames(vst_counts_B) %in% na_tray_samples]

pca_B_nona <- prcomp(t(vst_counts_B_nona_tray))
plot_pca_B_nona <- stats:::summary.prcomp(pca_B_nona)$importance[2, ]
barplot(plot_pca_B_nona, ylab = "Proportion of variance",
        ylim = c(0, 0.25), xlim = c(0, 21.7), las = 2)


pdf(file = "results/PCA_B_treat_nona.pdf", width = 10, height = 8)
ggbiplot(pca_B_nona, obs.scale = 1, var.scale = 1,
         groups = data_info_B_nona$treatment, ellipse = FALSE, var.axes = FALSE,
         circle = FALSE) + scale_shape_discrete(name = '') +
  geom_point((aes(color = data_info_B_nona$treatment)), size = 2) +
  theme_bw() + scale_color_manual(values = brewer.pal(4, "Set2"),
                                  name = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

pdf(file = "results/PCA_B_tray_nona.pdf", width = 10, height = 8)
ggbiplot(pca_B_nona, obs.scale = 1, var.scale = 1,
         groups = data_info_B_nona$tray_id, ellipse = FALSE, var.axes = FALSE,
         circle = FALSE) + scale_shape_discrete(name = '') +
  geom_point((aes(color = data_info_B_nona$tray_id)), size = 2) +
  theme_bw() + scale_color_discrete(name = "Tray ID") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

#Get expression distance between groups and among groups-------------------------------------------------
#Get metadata
metadataA <- data_info_A[, c("treatment", "tray_id")]
colnames(metadataA) <- c("Treatment", "Tray ID")

metadataB <- data_info_B[, c("treatment", "tray_id")]
colnames(metadataB) <- c("Treatment", "Tray ID")

sampleDists_A <- dist(t(assay(vsd_A)))
sampleDistMatrix_A <- as.matrix(sampleDists_A)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file = "results/dist_heatmapA.pdf", width = 10.5, height = 7.5)
pheatmap(sampleDistMatrix_A, annotation_col = metadataA,
         clustering_distance_rows = sampleDists_A,
         clustering_distance_cols = sampleDists_A,
         treeheight_row = 0, legend = FALSE,
         annotation_legend = TRUE, annotation_names_col = TRUE,
         col=colors)
dev.off()

sampleDists_B <- dist(t(assay(vsd_B)))
sampleDistMatrix_B <- as.matrix(sampleDists_B)
pdf(file = "results/dist_heatmapB.pdf", width = 10.5, height = 7.5)
pheatmap(sampleDistMatrix_B, annotation_col = metadataB,
         clustering_distance_rows = sampleDists_B,
         clustering_distance_cols = sampleDists_B,
         treeheight_row = 0, legend = FALSE,
         annotation_legend = TRUE, annotation_names_col = TRUE,
         col=colors)
dev.off()

#Distances in experiment A-------------------------------------------------------------------------------------------------
#Per sample, calculate its median distance to samples in other groups
expA_dists <- dist(t(assay(vsd_A)))
expA_dists_mat <- as.matrix(expA_dists)
all_median_dist_A <- matrix(NA, ncol = 4)

for (this_sample in colnames(expA_dists_mat)){
  dists_mat_sample <- expA_dists_mat[, this_sample]
  a1_samples <- data_info_A$sample_name[grepl("B1", data_info_A$treatment_raw)]
  a1_dist <- median(dists_mat_sample[a1_samples])
  a2_samples <- data_info_A$sample_name[grepl("A2", data_info_A$treatment_raw)]
  a2_dist <- median(dists_mat_sample[a2_samples])
  a3_samples <- data_info_A$sample_name[grepl("A3", data_info_A$treatment_raw)]
  a3_dist <- median(dists_mat_sample[a3_samples])
  all_median_dist_A <- rbind(all_median_dist_A, c(this_sample, a1_dist, a2_dist, a3_dist))
}

#Clean and add groups per sample
all_median_dist_A <- as.data.frame(all_median_dist_A)
all_median_dist_A <- all_median_dist_A[complete.cases(all_median_dist_A), ]
colnames(all_median_dist_A) <- c("sample_name", "a1", "a2", "a3")

all_median_dist_A <- left_join(all_median_dist_A, data_info_A)
all_median_dist_A$tray_id <- NULL
all_median_dist_A$merged <- NULL
all_median_dist_A$treatment_raw <- NULL

all_median_dist_A <- melt(all_median_dist_A, id.vars = c("sample_name", "treatment"))
all_median_dist_A$value <- as.numeric(all_median_dist_A$value)

#Plot the thing
pdf(file = "results/test_intragroupA.pdf", width = 10, height = 8)
ggplot(data = all_median_dist_A, aes(x = variable, y = value, color = treatment)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_dodge(0.8)) +
  xlab("Distance to samples in") + ylab("Median distance") + theme_bw() +
  facet_grid(. ~ variable, scales = "free") +
  scale_color_manual(values = brewer.pal(3, "Set1"),
                     name = "Treatment") +
  scale_x_discrete(labels = c("a1" = "NMQ", "a2" = "SFQ 3 dpf", "a3" = "GFQ 3 dpf")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))
dev.off()
#For B--------------
expB_dists <- dist(t(assay(vsd_B)))
expB_dists_mat <- as.matrix(expB_dists)
all_median_dist_b <- matrix(NA, ncol = 5)

for (this_sample in colnames(expB_dists_mat)){
  dists_mat_sample <- expB_dists_mat[, this_sample]
  b1_samples <- data_info_B$sample_name[grepl("B1", data_info_B$treatment_raw)]
  b1_dist <- median(dists_mat_sample[b1_samples])
  b2_samples <- data_info_B$sample_name[grepl("B2", data_info_B$treatment_raw)]
  b2_dist <- median(dists_mat_sample[b2_samples])
  b3_samples <- data_info_B$sample_name[grepl("B3", data_info_B$treatment_raw)]
  b3_dist <- median(dists_mat_sample[b3_samples])
  b4_samples <- data_info_B$sample_name[grepl("B4", data_info_B$treatment_raw)]
  b4_dist <- median(dists_mat_sample[b4_samples])
  all_median_dist_b <- rbind(all_median_dist_b, c(this_sample, b1_dist, b2_dist,
                                                  b3_dist, b4_dist))
}

#Clean and add groups per sample
all_median_dist_b <- as.data.frame(all_median_dist_b)
all_median_dist_b <- all_median_dist_b[complete.cases(all_median_dist_b), ]
colnames(all_median_dist_b) <- c("sample_name", "b1", "b2", "b3", "b4")

all_median_dist_b <- left_join(all_median_dist_b, data_info_B)
all_median_dist_b$tray_id <- NULL
all_median_dist_b$merged <- NULL
all_median_dist_b$treatment_raw <- NULL

all_median_dist_b <- melt(all_median_dist_b, id.vars = c("sample_name", "treatment"))
all_median_dist_b$value <- as.numeric(all_median_dist_b$value)

#Plot the thing
pdf(file = "results/test_intragroupB.pdf", width = 10, height = 8)
ggplot(data = all_median_dist_b, aes(x = variable, y = value, color = treatment)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_dodge(0.8)) +
  xlab("Distance to samples in") + ylab("Median distance") + theme_bw() +
  facet_grid(. ~ variable, scales = "free") +
  scale_color_manual(values = brewer.pal(4, "Set2"),
                     name = "Treatment") +
  scale_x_discrete(labels = c("b1" = "NMQ", "b2" = "SFQ 25dpf",
                              "b3" = "GFQsmall\n25dpf", "b4" = "GFQlarge\n25dpf")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))
dev.off()
