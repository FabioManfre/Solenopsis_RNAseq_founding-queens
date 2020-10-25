setwd("/Users/fabiomanfredini/Dropbox/SOLENOPSIS_RNAseq/ANALYSES/Carlos_NEW/NETWORK_analyses_NCBI/whole_geneset")

library(WGCNA)
options(stringsAsFactors = FALSE);
allowWGCNAThreads()


sinv.reads=read.table("raw_reads_whole-set_module-trait.txt",h=T, row.names=1)

head(sinv.reads)

     NMQ_1 P03_1 H03_1 P03_2 H03_2 P03_3 NMQ_2  NMQ_3 NMQ_4 H03_3 P03_4 P03_5 NMQ_5 H03_4 NMQ_6
ATP6  8956 10104 12536 11995 12360 10491 13065  17810 11856  7209  9180 12420  8666  8243 15161
ATP8  1122  1783  2260  2439  2736  2085  2079   4101  2471  1362  1677  2779  1875  2000  2612
COX1 77406 61457 97869 90305 77261 74774 90489 107396 78534 55534 73492 69636 63306 66549 96687
COX2  6948  9632 14261 13216  8710  9239 16010  12464  8779  6044 10323  9990 10118 12148 15807
COX3 10411 22038 22716 20912 22363 21976 25997  34063 23533 15102 15504 20786 15574 37241 21571
CYTB  9969 10175 17181 16429 15000 11988 16206  18371 11467 10015 14211 10884 10521 11196  9756

      P03_6 H03_5  H03_6 L25_1 L25_2 H25_1 H25_2 H25_3 H25_4 H25_5 S25_1 S25_2 L25_3 S25_3 L25_4
ATP6  23498  7439  72445  5640 12883 13402  7243  7182  5941 10218  8532 10081  9897  7127 13193
ATP8   5677  1876  18934   988  2583  2754  1456  1139  1438  2262  2589  3081  2153  1498  3789
COX1 135823 53198 387971 49753 68247 85063 57828 73299 41181 70193 50743 59854 76154 54245 75900
COX2  19276  9532  55413  7433  8773 14957  9015  8062  4521 10081  5778  6920 10991  5884  9158
COX3  39854 28926 122812 13053 18544 27543 13069 13598 11724 17168 15002 16958 17969 15610 18262
CYTB  23580 10985  79163 11180  9426 22790 11299 11780  7052 15863 10623 12767 16778 10612 11753

     S25_4 S25_5 L25_5
ATP6  6041 11954  2298
ATP8  1600  3149   491
COX1 46280 60223 38596
COX2  7810  8791  6326
COX3 22235 19011 21233
CYTB  9885 15651  5811

sinv.med=apply(sinv.reads,1,median)

head(sinv.med)

 ATP6  ATP8  COX1  COX2  COX3  CYTB 
10104  2153 70193  9239 20786 11467 

sinv.reads.tmp=sinv.reads[sinv.med>1,] 

sinv.reads=sinv.reads.tmp

datExpr.sinv= as.data.frame(t(log2(0.5+sinv.reads)));

head(datExpr.sinv)

      LOC105208173 LOC105208174 LOC105208176 LOC105208177 LOC105208178 LOC105208179 LOC105208181
NMQ_1    0.5849625     3.087463     5.189825     7.262095     4.977280     5.894818     2.700440
P03_1    2.1699250     4.285402     6.219169     7.909893     5.894818     6.238405     5.228822
H03_1    0.5849625     4.977280     6.219169     8.315150     6.467606     5.870365     5.022368
P03_2    1.3219281     4.209453     5.189825     6.982994     5.768184     5.266787     4.285402
H03_2    0.6770797     5.442943     5.942515     8.567956     6.562242     5.794416     5.228819
P03_3   -1.0000000     4.727920     6.011227     7.810572     5.988685     5.870365     5.228819

      LOC105208182 LOC105208183 LOC105208184 LOC105208185 LOC105208186 LOC105208187 LOC105208188
NMQ_1     7.224002     6.515700     7.071462     3.857981     9.565102     6.276124     7.154818
P03_1     7.661778     6.159871     7.344296     4.285402    10.025832     5.599913     7.353147
H03_1     8.147205     6.451211     7.873444     4.672425    10.484319     6.098032     7.891784
P03_2     7.299208     5.870365     7.317413     3.857981     9.439831     5.149747     6.948367
H03_2     8.216746     6.546894     8.350939     4.491853    10.576012     5.768184     7.968667
P03_3     7.558421     6.179909     7.764872     4.491853    10.018896     6.098032     7.535275



names(datExpr.sinv) = rownames(sinv.reads);



gsg.sinv= goodSamplesGenes(datExpr.sinv, verbose = 3);

 Flagging genes and samples with too many missing values...
  ..step 1

datExpr.sinv=datExpr.sinv[,gsg.sinv$goodGenes]

powers = c(c(1:10), seq(from = 12, to=30, by=2))

head (datExpr.sinv)

      LOC105208142 LOC105208143 LOC105208144 LOC105208145 LOC105208146 LOC105208147 LOC105208148
NMQ_1     4.426265     7.661778     8.564149     5.918863     9.591522     7.710807     5.022368
P03_1     5.022368     7.879583     9.251482     6.348728    10.025832     8.407265     3.643856
H03_1     4.614710     8.489848     9.365229     6.033423    10.085472     8.236014     3.247928
P03_2     4.357552     7.836050     8.911392     6.199672     9.265615     7.611025     4.209453
H03_2     5.507795     8.187352     9.696098     6.098032    10.297490     8.477764     3.857981
P03_3     5.189825     8.136991     9.360847     6.651052     9.941781     8.231221     4.209453

      LOC105208149 LOC105208150 LOC105208151 LOC105208152 LOC105208153 LOC105208154 LOC105208155
NMQ_1     6.330915     3.857981     7.422068     2.906891     9.861861     6.707359     5.741467
P03_1     6.366322     4.781360     8.407268     3.523562     9.655533     7.353147     4.977280
H03_1     6.864186     4.285402     8.292323     4.129283    10.152918     7.124121     5.686501
P03_2     5.686498     3.857981     7.945444     2.906891     9.195986     6.467606     5.022368
H03_2     6.960001     4.727920     8.445015     4.727920    10.358102     7.174926     5.442943
P03_3     6.761550     5.189825     8.453271     3.392317     9.633903     7.005625     5.409391


2.b Step-by-step network construction and module detection
2.b.1 Choosing the soft-thresholding power: analysis of network topology


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr.sinv, powerVector = powers, verbose = 5)

pickSoftThreshold: will use block size 3490.
 pickSoftThreshold: calculating connectivity for given powers...
   ..working on genes 1 through 3490 of 12817
   ..working on genes 3491 through 6980 of 12817
   ..working on genes 6981 through 10470 of 12817
   ..working on genes 10471 through 12817 of 12817
   Power SFT.R.sq    slope truncated.R.sq mean.k. median.k. max.k.
1      1 0.985000  2.48000          0.982    7820    8480.0   9940
2      2 0.923000  0.90900          0.903    5390    6010.0   8140
3      3 0.525000  0.30500          0.460    3970    4430.0   6870
4      4 0.000199 -0.00483         -0.212    3050    3350.0   5940
5      5 0.198000 -0.18300         -0.027    2420    2580.0   5200
6      6 0.385000 -0.30600          0.218    1960    2010.0   4610
7      7 0.496000 -0.40500          0.389    1610    1590.0   4130
8      8 0.558000 -0.48700          0.504    1350    1270.0   3720
9      9 0.605000 -0.55500          0.595    1140    1030.0   3370
10    10 0.626000 -0.62000          0.650     974     831.0   3070
11    12 0.651000 -0.74200          0.722     727     556.0   2590
12    14 0.677000 -0.83400          0.786     555     382.0   2210
13    16 0.700000 -0.91100          0.829     433     266.0   1900
14    18 0.723000 -0.97500          0.864     343     188.0   1650
15    20 0.741000 -1.04000          0.890     275     135.0   1450
16    22 0.750000 -1.09000          0.903     224      98.3   1270
17    24 0.757000 -1.15000          0.915     184      72.2   1130
18    26 0.772000 -1.19000          0.928     152      53.6   1000
19    28 0.779000 -1.23000          0.935     127      40.1    895
20    30 0.788000 -1.26000          0.942     107      30.2    803



# Plot the results:

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");


# this line corresponds to using an R^2 cut-off of h

abline(h=0.90,col="red")
abline(h=0.80,col="red")
abline(h=0.70,col="red")
abline(h=0.60,col="red")


# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

abline(h=1000,col="blue")
abline(h=1500,col="blue")
abline(h=500,col="blue")


###I choose 12 as power because it's the lowest power that doesn't go below 500 connectivity

BETA = 12

net.sinv12 = blockwiseModules(datExpr.sinv, power = BETA, minModuleSize = 10, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = F, verbose = 3)

Calculating module eigengenes block-wise from all genes
   Flagging genes and samples with too many missing values...
    ..step 1
 ....pre-clustering genes to determine blocks..
   Projective K-means:
   ..k-means clustering..
   ..merging smaller clusters...
Block sizes:
gBlocks
   1    2    3 
4999 4991 2827 
 ..Working on block 1 .
    TOM calculation: adjacency..
adjacency: replaceMissing: 0
    ..will use 8 parallel threads.
     Fraction of slow calculations: 0.000000
    ..connectivity..
    ..matrix multiplication..
    ..normalization..
    ..done.
 ....clustering..
 ....detecting modules..
 ....calculating module eigengenes..
 ....checking modules for statistical meaningfulness..
 ..Working on block 2 .
    TOM calculation: adjacency..
adjacency: replaceMissing: 0
    ..will use 8 parallel threads.
     Fraction of slow calculations: 0.000000
    ..connectivity..
    ..matrix multiplication..
    ..normalization..
    ..done.
 ....clustering..
 ....detecting modules..
 ....calculating module eigengenes..
 ....checking modules for statistical meaningfulness..
 ..Working on block 3 .
    TOM calculation: adjacency..
adjacency: replaceMissing: 0
    ..will use 8 parallel threads.
     Fraction of slow calculations: 0.000000
    ..connectivity..
    ..matrix multiplication..
    ..normalization..
    ..done.
 ....clustering..
 ....detecting modules..
 ....calculating module eigengenes..
 ....checking modules for statistical meaningfulness..
 ..merging modules that are too close..
     mergeCloseModules: Merging modules whose distance is less than 0.25
       Calculating new MEs...

 

sizeGrWindow(12, 9)
mergedColors = labels2colors(net.sinv12$colors)
plotDendroAndColors(net.sinv12$dendrograms[[1]], mergedColors[net.sinv12$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net.sinv12$colors;
moduleColors = labels2colors(net.sinv12$colors)
table(moduleColors)


moduleColors
    black      blue     brown     green      grey   magenta      pink    purple       red turquoise 
       20       110        26        22       435        15        16        15        20     12114 
       
   yellow 
       24 

############################################################################################
############################################################################################


TRAIT-MODULE ASSOCIATION STUDY

#Start with defining the design for the experiment

ReadCounts = read.table("raw_reads_whole-set_module-trait.txt", header = T, row.names = 1)

#Select the different treatment types by col names (in this case, the treatment
#is indicated by the digits 1-2 in the col name. The unique function then decreases
#the treatments even further to only those that are unique, in this example, there 
#were several samples for each treatment. Then changes the variable into a factor
#so that it can used by ANOVA

treatments = substr(colnames(ReadCounts), 1, 3)
treats = unique(treatments)
treatments = factor(treatments)
treatments 

 [1] NMQ P03 H03 P03 H03 P03 NMQ NMQ NMQ H03 P03 P03 NMQ H03 NMQ P03 H03 H03 L25 L25 H25 H25 H25 H25
 
[25] H25 S25 S25 L25 S25 L25 S25 S25 L25

Levels: H03 H25 L25 NMQ P03 S25

library(edgeR)
library(ggplot2)
library(RMySQL)
library(edgeR)
library(DESeq2)
library(GOstats)
library(GSEABase)
library(gplots)
library(gmodels)
library(gtools)
library(pgirmess)
library(class)
library(gridExtra)
library(reshape)
library(WGCNA)
library(Biobase)

#creates a normalization factor using the TMM method. 

normFact = calcNormFactors(ReadCounts, method = "TMM")
normFact

 [1] 1.1446879 0.9733246 0.9827481 0.8979541 1.0406718 0.9325592 1.1090251 0.9976741 1.0965625
[10] 0.9499181 1.0085856 0.9282652 1.1576250 0.8606962 1.0853980 0.9269470 1.0198586 1.0851393
[19] 0.9068517 1.0592514 0.8759787 0.9372158 1.0421639 1.0322595 1.1139614 0.8560495 0.9051030
[28] 0.9366407 1.1370559 1.1042829 0.9423225 0.9540250 1.1308167

#We create a DGEList object that holds the data, treatments and adjusted library sizes.

treatments = substr(colnames(ReadCounts), 1, 3)
d = DGEList(counts = ReadCounts, group = treatments, genes = rownames(ReadCounts), lib.size = colSums(ReadCounts) * normFact)
d

An object of class "DGEList"
$counts
     NMQ_1 P03_1 H03_1 P03_2 H03_2 P03_3 NMQ_2  NMQ_3 NMQ_4 H03_3 P03_4 P03_5 NMQ_5 H03_4 NMQ_6
ATP6  8956 10104 12536 11995 12360 10491 13065  17810 11856  7209  9180 12420  8666  8243 15161
ATP8  1122  1783  2260  2439  2736  2085  2079   4101  2471  1362  1677  2779  1875  2000  2612
COX1 77406 61457 97869 90305 77261 74774 90489 107396 78534 55534 73492 69636 63306 66549 96687
COX2  6948  9632 14261 13216  8710  9239 16010  12464  8779  6044 10323  9990 10118 12148 15807
COX3 10411 22038 22716 20912 22363 21976 25997  34063 23533 15102 15504 20786 15574 37241 21571

      P03_6 H03_5  H03_6 L25_1 L25_2 H25_1 H25_2 H25_3 H25_4 H25_5 S25_1 S25_2 L25_3 S25_3 L25_4
ATP6  23498  7439  72445  5640 12883 13402  7243  7182  5941 10218  8532 10081  9897  7127 13193
ATP8   5677  1876  18934   988  2583  2754  1456  1139  1438  2262  2589  3081  2153  1498  3789
COX1 135823 53198 387971 49753 68247 85063 57828 73299 41181 70193 50743 59854 76154 54245 75900
COX2  19276  9532  55413  7433  8773 14957  9015  8062  4521 10081  5778  6920 10991  5884  9158
COX3  39854 28926 122812 13053 18544 27543 13069 13598 11724 17168 15002 16958 17969 15610 18262

     S25_4 S25_5 L25_5
ATP6  6041 11954  2298
ATP8  1600  3149   491
COX1 46280 60223 38596
COX2  7810  8791  6326
COX3 22235 19011 21233

14461 more rows ...

$samples
      group lib.size norm.factors
NMQ_1   NMQ  4817063            1
P03_1   P03  5825867            1
H03_1   H03  6255457            1
P03_2   P03  4024500            1
H03_2   H03  7261957            1
28 more rows ...

$genes
  genes
1  ATP6
2  ATP8
3  COX1
4  COX2
5  COX3
14461 more rows ...


d$samples

      group lib.size norm.factors
NMQ_1   NMQ  4817063            1
P03_1   P03  5825867            1
H03_1   H03  6255457            1
P03_2   P03  4024500            1
H03_2   H03  7261957            1
P03_3   P03  5678273            1
NMQ_2   NMQ  6761558            1
NMQ_3   NMQ  7777177            1
NMQ_4   NMQ  7069651            1
H03_3   H03  5860432            1
P03_4   P03  5620128            1
P03_5   P03  4697766            1
NMQ_5   NMQ  5639013            1
H03_4   H03  6193382            1
NMQ_6   NMQ  5471461            1
P03_6   P03  8091996            1
H03_5   H03  6650916            1
H03_6   H03 48548547            1
L25_1   L25  4279824            1
L25_2   L25  6182327            1
H25_1   H25  7269743            1
H25_2   H25  5823905            1
H25_3   H25  6075334            1
H25_4   H25  6023240            1
H25_5   H25  7894884            1
S25_1   S25  4997606            1
S25_2   S25  4816628            1
L25_3   L25  5079008            1
S25_3   S25  6297283            1
L25_4   L25  6886559            1
S25_4   S25  6062290            1
S25_5   S25  5815788            1
L25_5   L25  5936472            1

levels(d$samples$group)

[1] "H03" "H25" "L25" "NMQ" "P03" "S25"


#create design matrix, and label rows and columns according to our contrasts

design <- model.matrix(~0+group, data=d$samples)
colnames(design) <- levels(d$samples$group)
design

      H03 H25 L25 NMQ P03 S25
NMQ_1   0   0   0   1   0   0
P03_1   0   0   0   0   1   0
H03_1   1   0   0   0   0   0
P03_2   0   0   0   0   1   0
H03_2   1   0   0   0   0   0
P03_3   0   0   0   0   1   0
NMQ_2   0   0   0   1   0   0
NMQ_3   0   0   0   1   0   0
NMQ_4   0   0   0   1   0   0
H03_3   1   0   0   0   0   0
P03_4   0   0   0   0   1   0
P03_5   0   0   0   0   1   0
NMQ_5   0   0   0   1   0   0
H03_4   1   0   0   0   0   0
NMQ_6   0   0   0   1   0   0
P03_6   0   0   0   0   1   0
H03_5   1   0   0   0   0   0
H03_6   1   0   0   0   0   0
L25_1   0   0   1   0   0   0
L25_2   0   0   1   0   0   0
H25_1   0   1   0   0   0   0
H25_2   0   1   0   0   0   0
H25_3   0   1   0   0   0   0
H25_4   0   1   0   0   0   0
H25_5   0   1   0   0   0   0
S25_1   0   0   0   0   0   1
S25_2   0   0   0   0   0   1
L25_3   0   0   1   0   0   0
S25_3   0   0   0   0   0   1
L25_4   0   0   1   0   0   0
S25_4   0   0   0   0   0   1
S25_5   0   0   0   0   0   1
L25_5   0   0   1   0   0   0
attr(,"assign")
[1] 1 1 1 1 1 1
attr(,"contrasts")
attr(,"contrasts")$group
[1] "contr.treatment"


***********************************
***********************************

# Recalculate MEs with color labels

invisible(MEs0 <- moduleEigengenes(datExpr.sinv, moduleColors)$eigengenes)
MEs = orderMEs(MEs0)
MEs<-MEs[,order(names(MEs))]

##ALL QUEEN PHENOTYPES SEPARATE
###This is a quite convoluted way to change the order of the queen groups in the picture...I couldn't manage to do it in a different way :-()

traits_sorted = cbind(NMQs=design[,"NMQ"], single03=design[,"H03"], group03=design[,"P03"], single25=design[,"H25"], small25=design[,"S25"], large25=design[,"L25"])
traits_sorted


moduleTraitCor = cor(MEs, traits_sorted, use = "p")
moduleTraitPvalue = p.adjust(corPvalueStudent(moduleTraitCor, nSamples=nrow(MEs)),method="fdr")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = colnames(moduleTraitCor),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = .5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))


mycols <- colorpanel(n=19,low="white",high="darkgreen")
heatmap.2(moduleTraitCor,key=FALSE,trace="none",col=mycols,Colv=NULL,labCol=colnames(traits_sorted),
density.info="none", hclustfun=function(c){hclust(c, method="complete")},margins = c(10, 4), 
RowSideColors=gsub("ME","",names(MEs)),labRow=table(moduleColors))


##Plotting module-all-queens correlation coefficient over time

moduleTraitCor_stack <- melt(moduleTraitCor)
colnames(moduleTraitCor_stack) <- c("module","stage","corr")
moduleTraitCor_stack$stage <- factor(moduleTraitCor_stack$stage,unique(moduleTraitCor_stack[,2]))
ggplot(subset(moduleTraitCor_stack, stage !="single" & stage != "single_5"),aes(x=stage,y=corr,
color=module,group=module)) + geom_line()+theme_bw()



##Plotting module-group correlation coefficient over time

moduleTraitCor_stack <- melt(moduleTraitCor)
colnames(moduleTraitCor_stack) <- c("module","stage","corr")
moduleTraitCor_stack$stage <- factor(moduleTraitCor_stack$stage,unique(moduleTraitCor_stack[,2]))
ggplot(subset(moduleTraitCor_stack, stage !="single03" & stage !="single25"),aes(x=stage,y=corr,
color=module,group=module)) + geom_line()+theme_bw()


##Plotting module-single correlation coefficient over time

moduleTraitCor_stack <- melt(moduleTraitCor)
colnames(moduleTraitCor_stack) <- c("module","stage","corr")
moduleTraitCor_stack$stage <- factor(moduleTraitCor_stack$stage,unique(moduleTraitCor_stack[,2]))
ggplot(subset(moduleTraitCor_stack, stage !="group03" & stage != "small25" & stage != "large25"),aes(x=stage,y=corr,
color=module,group=module)) + geom_line()+theme_bw()


##Plotting module-3_days correlation coefficient

moduleTraitCor_stack <- melt(moduleTraitCor)
colnames(moduleTraitCor_stack) <- c("module","stage","corr")
moduleTraitCor_stack$stage <- factor(moduleTraitCor_stack$stage,unique(moduleTraitCor_stack[,2]))
ggplot(subset(moduleTraitCor_stack, stage !="single25" & stage !="small25" & stage !="large25"),aes(x=stage,y=corr,
color=module,group=module)) + geom_line()+theme_bw()


##Plotting module-25_days correlation coefficient

moduleTraitCor_stack <- melt(moduleTraitCor)
colnames(moduleTraitCor_stack) <- c("module","stage","corr")
moduleTraitCor_stack$stage <- factor(moduleTraitCor_stack$stage,unique(moduleTraitCor_stack[,2]))
ggplot(subset(moduleTraitCor_stack, stage !="single03" & stage !="group03"),aes(x=stage,y=corr,
color=module,group=module)) + geom_line()+theme_bw()

>>these plots are a sort of K-means clustering: for each module, they show how its correlation with the different queen phenotypes goes up or down, hence how the patterns of genes in each module change over time/modality of colony founding!

************************************

##COMBINED PHENOTYPES
group vs. single
3 days vs. 25 days
founding vs. non-founding

traits_grouped = cbind(group=design[,"P03"]+design[,"S25"]+design[,"L25"], single=design[,"H03"]+design[,"H25"], early=design[,"H03"]+design[,"P03"], late=design[,"H25"]+design[,"S25"]+design[,"L25"], foundress=design[,"H03"]+design[,"P03"]+design[,"H25"]+design[,"S25"]+design[,"L25"], non_foundress=design[,"NMQ"])

moduleTraitCor = cor(MEs, traits_grouped, use = "p")
moduleTraitPvalue = p.adjust(corPvalueStudent(moduleTraitCor, nSamples=nrow(MEs)),method="fdr")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = colnames(moduleTraitCor),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = .5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))



mycols <- colorpanel(n=19,low="white",high="darkgreen")

heatmap.2(moduleTraitCor,key=FALSE,trace="none",col=mycols,Colv=FALSE,dendrogram="row", labCol=colnames(traits_grouped),
density.info="none", hclustfun=function(c){hclust(c, method="complete")},margins = c(10, 4),
RowSideColors=gsub("ME","",names(MEs)),labRow=table(moduleColors))

******************************************************************************************

6 Exporting network data to network visualization software
6.a Exporting to VisANT

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr.sinv, power = 12);

##STANDARD PROTOCOL

setwd("/Users/fabiomanfredini/Dropbox/SOLENOPSIS_RNAseq/ANALYSES/Carlos_NEW/NETWORK_analyses_NCBI/whole_geneset/modules")

###GREEN
# Select module
module = "green";


# Select module probes
probes = names(datExpr.sinv)
inModule = (moduleColors==module);
modProbes = probes[inModule];


# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)


# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

###PURPLE
# Select module
module = "purple";


# Select module probes
probes = names(datExpr.sinv)
inModule = (moduleColors==module);
modProbes = probes[inModule];


# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)


# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTinput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

###BLUE
# Select module
module = "blue";


# Select module probes
probes = names(datExpr.sinv)
inModule = (moduleColors==module);
modProbes = probes[inModule];


# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)


# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

###BLACK
# Select module
module = "black";


# Select module probes
probes = names(datExpr.sinv)
inModule = (moduleColors==module);
modProbes = probes[inModule];


# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)


# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

###PINK
# Select module
module = "pink";


# Select module probes
probes = names(datExpr.sinv)
inModule = (moduleColors==module);
modProbes = probes[inModule];


# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)


# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

###MAGENTA
# Select module
module = "magenta";


# Select module probes
probes = names(datExpr.sinv)
inModule = (moduleColors==module);
modProbes = probes[inModule];


# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)


# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0)

