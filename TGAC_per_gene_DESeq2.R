# Summarise counts per gene (rather than transcript)
# then use for DESeq2
#
# Philippa Borrill
# 15/08/2016

# Kallisto's default output is per transcript but it is probably more informative
# to carry out DE analysis using genes rather than transcripts

# found a package called "tximportData" which should sum the transcript counts per gene for kallisto output
# details here: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

#setwd("//nbi//group-data//ifs//NBI//Research-Groups//Cristobal-Uauy//PB_AFLF//control_timecourse//TGAC_kallisto_analysis//kallisto_results_bootstrap//results/")
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results")

source("https://bioconductor.org/biocLite.R")
biocLite("tximportData")
install.packages("readr")
library(tximportData)
library(readr)

# make vector pointing to the kallisto results files
samples <- read.table("samples.txt", header=T)
samples

files <- file.path("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results", samples$sample, "abundance.tsv", fsep ="\\")
files
names(files) <- paste0(samples$sample)
files
all(file.exists(files))

# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.table("transcripts_to_genes.txt", header=T)
head(tx2gene)

library(tximport)
# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
names(txi)

head(txi$counts)
colnames(txi$counts)
# therefore txi$counts contains counts per gene for each sample


####### DESEQ2 analysis #########
# now use these aggregated counts for DESeq2

library(DESeq2)
timepoints <- c("FLB3", "FLB7", "FLB10", "FLB13", "FLB15", "FLB17", "FLB19", "FLB21", "FLB23","FLB26","G3", "G7", "G10", "G13", "G15", "G17", "G19", "G21", "G23","G26")
CondVector <- rep(timepoints,each=3)
CondVector

sampleTable <- data.frame(row.names=colnames(txi$counts), condition=as.factor(CondVector))
sampleTable

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

#make sure the FLB3 is used as the reference condition :
dds$condition <- relevel(dds$condition, "FLB3")

dim(dds)
# filter to remove genes with only 0 or 1 read
dds <- dds[ rowSums(counts(dds)) > 1, ]
dim(dds)

# saved a copy of dds just in case
dds_copy <- dds

# run DeSeq2
dds <- DESeq(dds)

# will need to have the bck_CDS1 set correctly to match old script which used this name for dds
#uses tutorial from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/ tutorial:
bck_CDS1 <- dds


#try out using http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/ tutorial:

#transform the raw counts so we can do clustering:
rld <- rlogTransformation(dds, blind=TRUE) ## can't do rld because too many genes therefore calc too slow - tutorial recommends to use vsd in this case

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

#plot graph to see effect of transformation #### WONT work 15th AUG 2016
pdf (file="DESeq2_VST_and_log2.pdf", height=5, width = 5)
par(mai=c(1,1,1,1))
px     <- counts(dds)[,1] / estimateSizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright", legend = c(expression("variance stabilizing transformation"), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

#plot graph to show comparison of transformation methods  #### WONT work 15th AUG 2016
library("vsn")
pdf (file="DESeq2_comparison_of_transformations.pdf", height=8, width = 8)
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,2.5))
meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))
meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))
dev.off()



#plot heatmap of samples

library("RColorBrewer")
library("gplots")


select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

pdf(file="DESeq2_heatmap1raw.pdf", height=8, width = 8)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
dev.off()


pdf(file="DESeq2_heatmap3varstab.pdf", height=8, width = 8)
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()






#try with vsd (rather than rld) to look at clustering 

distsVSD <- dist(t(assay(vsd)))
mat <- as.matrix(distsVSD)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))



pdf(file="deseq2_heatmaps_samplebysample_VSD.pdf", height =8, width =8)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13), cexRow=0.5, cexCol=0.5)
dev.off()


### separate out FLB and G
head(mat)
dim(mat)
FLB_mat<- mat[1:30,1:30]
dim(FLB_mat)

pdf(file="deseq2_heatmaps_samplebysample_vsd_FLB.pdf", height =8, width =8)
heatmap.2(FLB_mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

G_mat<- mat[31:60,31:60]
dim(G_mat)

pdf(file="deseq2_heatmaps_samplebysample_vsd_G.pdf", height =8, width =8)
heatmap.2(G_mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

#try PCA 
pdf(file="deseq2_pca_vsd.pdf", height =8, width =10)
print(plotPCA(vsd, intgroup=c("condition")))

dev.off()

# want separate pca for FLB and G. 

# 1st separate vsd
FLB_vsd <- vsd[,1:30]
pdf(file="deseq2_pca_vsd_FLB.pdf", height =8, width =10)
print(plotPCA(FLB_vsd, intgroup=c("condition")))
text()
dev.off()

# try plotting custom plot with ggplot
library(ggplot2)
pdf(file="deseq2_pca_vsd_FLB_with_names_full.pdf", height =8, width =10)
pca_FLB_data <- plotPCA(FLB_vsd, intgroup=c("condition"), returnData=TRUE)
ggplot(pca_FLB_data, aes(PC1, PC2, color=condition)) + geom_point(size=3) +
  geom_text(aes(label=row.names(pca_FLB_data)),hjust=-0.2, vjust=0.2)
dev.off()

pdf(file="deseq2_pca_vsd_FLB_with_names_cond.pdf", height =8, width =10)
pca_FLB_data <- plotPCA(FLB_vsd, intgroup=c("condition"), returnData=TRUE)
ggplot(pca_FLB_data, aes(PC1, PC2, color=condition)) + geom_point(size=3) +
  geom_text(aes(label=condition),hjust=-0.2, vjust=0.2)
dev.off()

# 1st separate vsd
G_vsd <- vsd[,31:60]
pdf(file="deseq2_pca_vsd_G.pdf", height =8, width =10)
print(plotPCA(G_vsd, intgroup=c("condition")))

dev.off()


# make PCA with colour coding
#data <- plotPCA(rld, intgroup=c("condition"),returnData=TRUE)
#percentVar <- round(100 * attr(data, "percentVar"))

#library("ggplot2")

#ggplot(data, aes(PC1, PC2, color=condition)) + geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
# ylab(paste0("PC2: ",percentVar[2],"% variance"))



##### got to here 9.2.16 ########

# want to save the results of DESeq (i.e. bckCDS_1)
# all timepoints + tissues
timepoints

FLB_timepoints <- timepoints[1:10]
FLB_timepoints
FLB_timepoints_no_FLB3 <- timepoints[2:10]

G_timepoints <- timepoints[11:20]
G_timepoints
G_timepoints_no_G3 <-timepoints[12:20]


# get data for all G timepoints vs G3
for (i in G_timepoints_no_G3) { 
  
  bck_res <- results(bckCDS_1,contrast=c("condition",i,"G3"))
  
  # sort results on padj
  ordered_res <- bck_res[order(bck_res$padj),]
  head(ordered_res)
  
  # output ordered_res to csv
  write.csv(ordered_res,file=paste(i,"_vs_G3_results.csv"))
  
  ordered_res_na.rm <- na.omit(ordered_res)
  
  assign(paste(i,"_vs_G3_results_na.rm",sep=""), ordered_res_na.rm[ordered_res_na.rm$padj<0.05,]) 
  
  print(paste(i,"_vs_G3_results_na.rm",sep=""))
  dim(assign(paste(i,"_vs_G3_results_na.rm",sep=""), ordered_res_na.rm[ordered_res_na.rm$padj<0.05,]) )
  
}


# get data for all FLB timepoints vs FLB3
for (i in FLB_timepoints_no_FLB3) { 
 
  bck_res <- results(bckCDS_1,contrast=c("condition",i,"FLB3"))
  
  # sort results on padj
  ordered_res <- bck_res[order(bck_res$padj),]
  head(ordered_res)
  
  # output ordered_res to csv
  write.csv(ordered_res,file=paste(i,"_vs_FLB3_results.csv"))
  
  ordered_res_na.rm <- na.omit(ordered_res)
  
  assign(paste(i,"_vs_FLB3_results_na.rm",sep=""), ordered_res_na.rm[ordered_res_na.rm$padj<0.05,]) 
  
  print(paste(i,"_vs_FLB3_results_na.rm",sep=""))
  dim(assign(paste(i,"_vs_FLB3_results_na.rm",sep=""), ordered_res_na.rm[ordered_res_na.rm$padj<0.05,]) )
  
}

# for FLB compare the timepoint to the previous one

for (k in 1:9) {
  
#  print(FLB_timepoints[i])
#  print(FLB_timepoints[(i+1)])
  
  j <- FLB_timepoints[k]
  i <- FLB_timepoints[(k+1)]

  
# get data for timepoint vs. previous timepoint
#i <- "G26"
#j <- "G23"

bck_res <- results(bckCDS_1,contrast=c("condition",i,j))

# sort results on padj
ordered_res <- bck_res[order(bck_res$padj),]
head(ordered_res)

# output ordered_res to csv
write.csv(ordered_res,file=paste(i,"_vs_",j,"_results.csv"))
ordered_res_na_rm <- na.omit(ordered_res)

print(paste(i,"_vs_",j))

# number of genes padj < 0.01
print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01,]))

# check number of genes DE padj<0.01 and 2x down reg
print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01 & ordered_res_na_rm$log2FoldChange<0.5,]))

# check number of genes DE padj<0.01 and 2x up reg
print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01 & ordered_res_na_rm$log2FoldChange>2,]))

# number of genes padj < 0.001
print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.001,]))

# check number of genes DE padj<0.01 and 2x down reg
print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.001 & ordered_res_na_rm$log2FoldChange<0.5,]))

# check number of genes DE padj<0.01 and 2x up reg
print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.001 & ordered_res_na_rm$log2FoldChange>2,]))
} 

# for G compare the timepoint to the previous one

for (k in 1:9) {
  
  #  print(G_timepoints[i])
  #  print(G_timepoints[(i+1)])
  
  j <- G_timepoints[k]
  i <- G_timepoints[(k+1)]
  
  
  bck_res <- results(bckCDS_1,contrast=c("condition",i,j))
  
  # sort results on padj
  ordered_res <- bck_res[order(bck_res$padj),]
  head(ordered_res)
  
  # output ordered_res to csv
  write.csv(ordered_res,file=paste(i,"_vs_",j,"_results.csv"))
  ordered_res_na_rm <- na.omit(ordered_res)
  
  print(paste(i,"_vs_",j))
  
  # number of genes padj < 0.01
  print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01,]))
  
  # check number of genes DE padj<0.01 and 2x down reg
  print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01 & ordered_res_na_rm$log2FoldChange<0.5,]))
  
  # check number of genes DE padj<0.01 and 2x up reg
  print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01 & ordered_res_na_rm$log2FoldChange>2,]))
  
  # number of genes padj < 0.001
  print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.001,]))
  
  # check number of genes DE padj<0.01 and 2x down reg
  print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.001 & ordered_res_na_rm$log2FoldChange<0.5,]))
  
  # check number of genes DE padj<0.01 and 2x up reg
  print(dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.001 & ordered_res_na_rm$log2FoldChange>2,]))
} 

#check number of genes DE
dim(FLB7_vs_FLB3_results_na.rm)
dim(FLB10_vs_FLB3_results_na.rm)
dim(FLB13_vs_FLB3_results_na.rm)
dim(FLB17_vs_FLB3_results_na.rm)
dim(FLB15_vs_FLB3_results_na.rm)
dim(FLB19_vs_FLB3_results_na.rm)
dim(FLB21_vs_FLB3_results_na.rm)
dim(FLB23_vs_FLB3_results_na.rm)
dim(FLB26_vs_FLB3_results_na.rm)

dim(G7_vs_G3_results_na.rm)
dim(G10_vs_G3_results_na.rm)
dim(G13_vs_G3_results_na.rm)
dim(G15_vs_G3_results_na.rm)
dim(G17_vs_G3_results_na.rm)
dim(G19_vs_G3_results_na.rm)
dim(G21_vs_G3_results_na.rm)
dim(G23_vs_G3_results_na.rm)
dim(G26_vs_G3_results_na.rm)


# check number of genes DE padj<0.01

dim(FLB7_vs_FLB3_results_na.rm[FLB7_vs_FLB3_results_na.rm$padj<0.01,])
dim(FLB10_vs_FLB3_results_na.rm[FLB10_vs_FLB3_results_na.rm$padj<0.01,])
dim(FLB13_vs_FLB3_results_na.rm[FLB13_vs_FLB3_results_na.rm$padj<0.01,])
dim(FLB15_vs_FLB3_results_na.rm[FLB15_vs_FLB3_results_na.rm$padj<0.01,])
dim(FLB17_vs_FLB3_results_na.rm[FLB17_vs_FLB3_results_na.rm$padj<0.01,])
dim(FLB19_vs_FLB3_results_na.rm[FLB19_vs_FLB3_results_na.rm$padj<0.01,])
dim(FLB21_vs_FLB3_results_na.rm[FLB21_vs_FLB3_results_na.rm$padj<0.01,])
dim(FLB23_vs_FLB3_results_na.rm[FLB23_vs_FLB3_results_na.rm$padj<0.01,])
dim(FLB26_vs_FLB3_results_na.rm[FLB26_vs_FLB3_results_na.rm$padj<0.01,])

dim(G7_vs_G3_results_na.rm[G7_vs_G3_results_na.rm$padj<0.01,])
dim(G10_vs_G3_results_na.rm[G10_vs_G3_results_na.rm$padj<0.01,])
dim(G13_vs_G3_results_na.rm[G13_vs_G3_results_na.rm$padj<0.01,])
dim(G15_vs_G3_results_na.rm[G15_vs_G3_results_na.rm$padj<0.01,])
dim(G17_vs_G3_results_na.rm[G17_vs_G3_results_na.rm$padj<0.01,])
dim(G19_vs_G3_results_na.rm[G19_vs_G3_results_na.rm$padj<0.01,])
dim(G21_vs_G3_results_na.rm[G21_vs_G3_results_na.rm$padj<0.01,])
dim(G23_vs_G3_results_na.rm[G23_vs_G3_results_na.rm$padj<0.01,])
dim(G26_vs_G3_results_na.rm[G26_vs_G3_results_na.rm$padj<0.01,])

# check number of genes DE padj<0.01 and 2x down reg
dim(FLB7_vs_FLB3_results_na.rm[FLB7_vs_FLB3_results_na.rm$padj<0.01 & FLB7_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB10_vs_FLB3_results_na.rm[FLB10_vs_FLB3_results_na.rm$padj<0.01 & FLB10_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB13_vs_FLB3_results_na.rm[FLB13_vs_FLB3_results_na.rm$padj<0.01 & FLB13_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB15_vs_FLB3_results_na.rm[FLB15_vs_FLB3_results_na.rm$padj<0.01 & FLB15_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB17_vs_FLB3_results_na.rm[FLB17_vs_FLB3_results_na.rm$padj<0.01 & FLB17_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB19_vs_FLB3_results_na.rm[FLB19_vs_FLB3_results_na.rm$padj<0.01 & FLB19_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB21_vs_FLB3_results_na.rm[FLB21_vs_FLB3_results_na.rm$padj<0.01 & FLB21_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB23_vs_FLB3_results_na.rm[FLB23_vs_FLB3_results_na.rm$padj<0.01 & FLB23_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB26_vs_FLB3_results_na.rm[FLB26_vs_FLB3_results_na.rm$padj<0.01 & FLB26_vs_FLB3_results_na.rm$log2FoldChange<0.5,])

dim(G7_vs_G3_results_na.rm[G7_vs_G3_results_na.rm$padj<0.01 & G7_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G10_vs_G3_results_na.rm[G10_vs_G3_results_na.rm$padj<0.01 & G10_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G13_vs_G3_results_na.rm[G13_vs_G3_results_na.rm$padj<0.01 & G13_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G15_vs_G3_results_na.rm[G15_vs_G3_results_na.rm$padj<0.01 & G15_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G17_vs_G3_results_na.rm[G17_vs_G3_results_na.rm$padj<0.01 & G17_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G19_vs_G3_results_na.rm[G19_vs_G3_results_na.rm$padj<0.01 & G19_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G21_vs_G3_results_na.rm[G21_vs_G3_results_na.rm$padj<0.01 & G21_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G23_vs_G3_results_na.rm[G23_vs_G3_results_na.rm$padj<0.01 & G23_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G26_vs_G3_results_na.rm[G26_vs_G3_results_na.rm$padj<0.01 & G26_vs_G3_results_na.rm$log2FoldChange<0.5,])

# check number of genes DE padj<0.01 and 2x up reg
dim(FLB7_vs_FLB3_results_na.rm[FLB7_vs_FLB3_results_na.rm$padj<0.01 & FLB7_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB10_vs_FLB3_results_na.rm[FLB10_vs_FLB3_results_na.rm$padj<0.01 & FLB10_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB13_vs_FLB3_results_na.rm[FLB13_vs_FLB3_results_na.rm$padj<0.01 & FLB13_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB15_vs_FLB3_results_na.rm[FLB15_vs_FLB3_results_na.rm$padj<0.01 & FLB15_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB17_vs_FLB3_results_na.rm[FLB17_vs_FLB3_results_na.rm$padj<0.01 & FLB17_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB19_vs_FLB3_results_na.rm[FLB19_vs_FLB3_results_na.rm$padj<0.01 & FLB19_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB21_vs_FLB3_results_na.rm[FLB21_vs_FLB3_results_na.rm$padj<0.01 & FLB21_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB23_vs_FLB3_results_na.rm[FLB23_vs_FLB3_results_na.rm$padj<0.01 & FLB23_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB26_vs_FLB3_results_na.rm[FLB26_vs_FLB3_results_na.rm$padj<0.01 & FLB26_vs_FLB3_results_na.rm$log2FoldChange>2,])

dim(G7_vs_G3_results_na.rm[G7_vs_G3_results_na.rm$padj<0.01 & G7_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G10_vs_G3_results_na.rm[G10_vs_G3_results_na.rm$padj<0.01 & G10_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G13_vs_G3_results_na.rm[G13_vs_G3_results_na.rm$padj<0.01 & G13_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G15_vs_G3_results_na.rm[G15_vs_G3_results_na.rm$padj<0.01 & G15_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G17_vs_G3_results_na.rm[G17_vs_G3_results_na.rm$padj<0.01 & G17_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G19_vs_G3_results_na.rm[G19_vs_G3_results_na.rm$padj<0.01 & G19_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G21_vs_G3_results_na.rm[G21_vs_G3_results_na.rm$padj<0.01 & G21_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G23_vs_G3_results_na.rm[G23_vs_G3_results_na.rm$padj<0.01 & G23_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G26_vs_G3_results_na.rm[G26_vs_G3_results_na.rm$padj<0.01 & G26_vs_G3_results_na.rm$log2FoldChange>2,])




# check number of genes DE padj<0.01 and 2x up or down reg
dim(FLB7_vs_FLB3_results_na.rm[FLB7_vs_FLB3_results_na.rm$padj<0.01 & (FLB7_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB7_vs_FLB3_results_na.rm$log2FoldChange>2),])
dim(FLB10_vs_FLB3_results_na.rm[FLB10_vs_FLB3_results_na.rm$padj<0.01 & (FLB10_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB10_vs_FLB3_results_na.rm$log2FoldChange>2),])
dim(FLB13_vs_FLB3_results_na.rm[FLB13_vs_FLB3_results_na.rm$padj<0.01 & (FLB13_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB13_vs_FLB3_results_na.rm$log2FoldChange>2),])
dim(FLB15_vs_FLB3_results_na.rm[FLB15_vs_FLB3_results_na.rm$padj<0.01 & (FLB15_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB15_vs_FLB3_results_na.rm$log2FoldChange>2),])
dim(FLB17_vs_FLB3_results_na.rm[FLB17_vs_FLB3_results_na.rm$padj<0.01 & (FLB17_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB17_vs_FLB3_results_na.rm$log2FoldChange>2),])
dim(FLB19_vs_FLB3_results_na.rm[FLB19_vs_FLB3_results_na.rm$padj<0.01 & (FLB19_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB19_vs_FLB3_results_na.rm$log2FoldChange>2),])
dim(FLB21_vs_FLB3_results_na.rm[FLB21_vs_FLB3_results_na.rm$padj<0.01 & (FLB21_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB21_vs_FLB3_results_na.rm$log2FoldChange>2),])
dim(FLB23_vs_FLB3_results_na.rm[FLB23_vs_FLB3_results_na.rm$padj<0.01 & (FLB23_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB23_vs_FLB3_results_na.rm$log2FoldChange>2),])
dim(FLB26_vs_FLB3_results_na.rm[FLB26_vs_FLB3_results_na.rm$padj<0.01 & (FLB26_vs_FLB3_results_na.rm$log2FoldChange<0.5|FLB26_vs_FLB3_results_na.rm$log2FoldChange>2),])

dim(G7_vs_G3_results_na.rm[G7_vs_G3_results_na.rm$padj<0.01 & (G7_vs_G3_results_na.rm$log2FoldChange<0.5|G7_vs_G3_results_na.rm$log2FoldChange>2),])
dim(G10_vs_G3_results_na.rm[G10_vs_G3_results_na.rm$padj<0.01 & (G10_vs_G3_results_na.rm$log2FoldChange<0.5|G10_vs_G3_results_na.rm$log2FoldChange>2),])
dim(G13_vs_G3_results_na.rm[G13_vs_G3_results_na.rm$padj<0.01 & (G13_vs_G3_results_na.rm$log2FoldChange<0.5|G13_vs_G3_results_na.rm$log2FoldChange>2),])
dim(G15_vs_G3_results_na.rm[G15_vs_G3_results_na.rm$padj<0.01 & (G15_vs_G3_results_na.rm$log2FoldChange<0.5|G15_vs_G3_results_na.rm$log2FoldChange>2),])
dim(G17_vs_G3_results_na.rm[G17_vs_G3_results_na.rm$padj<0.01 & (G17_vs_G3_results_na.rm$log2FoldChange<0.5|G17_vs_G3_results_na.rm$log2FoldChange>2),])
dim(G19_vs_G3_results_na.rm[G19_vs_G3_results_na.rm$padj<0.01 & (G19_vs_G3_results_na.rm$log2FoldChange<0.5|G19_vs_G3_results_na.rm$log2FoldChange>2),])
dim(G21_vs_G3_results_na.rm[G21_vs_G3_results_na.rm$padj<0.01 & (G21_vs_G3_results_na.rm$log2FoldChange<0.5|G21_vs_G3_results_na.rm$log2FoldChange>2),])
dim(G23_vs_G3_results_na.rm[G23_vs_G3_results_na.rm$padj<0.01 & (G23_vs_G3_results_na.rm$log2FoldChange<0.5|G23_vs_G3_results_na.rm$log2FoldChange>2),])
dim(G26_vs_G3_results_na.rm[G26_vs_G3_results_na.rm$padj<0.01 & (G26_vs_G3_results_na.rm$log2FoldChange<0.5|G26_vs_G3_results_na.rm$log2FoldChange>2),])


# check number of genes DE padj<0.001

dim(FLB7_vs_FLB3_results_na.rm[FLB7_vs_FLB3_results_na.rm$padj<0.001,])
dim(FLB10_vs_FLB3_results_na.rm[FLB10_vs_FLB3_results_na.rm$padj<0.001,])
dim(FLB13_vs_FLB3_results_na.rm[FLB13_vs_FLB3_results_na.rm$padj<0.001,])
dim(FLB15_vs_FLB3_results_na.rm[FLB15_vs_FLB3_results_na.rm$padj<0.001,])
dim(FLB17_vs_FLB3_results_na.rm[FLB17_vs_FLB3_results_na.rm$padj<0.001,])
dim(FLB19_vs_FLB3_results_na.rm[FLB19_vs_FLB3_results_na.rm$padj<0.001,])
dim(FLB21_vs_FLB3_results_na.rm[FLB21_vs_FLB3_results_na.rm$padj<0.001,])
dim(FLB23_vs_FLB3_results_na.rm[FLB23_vs_FLB3_results_na.rm$padj<0.001,])
dim(FLB26_vs_FLB3_results_na.rm[FLB26_vs_FLB3_results_na.rm$padj<0.001,])

dim(G7_vs_G3_results_na.rm[G7_vs_G3_results_na.rm$padj<0.001,])
dim(G10_vs_G3_results_na.rm[G10_vs_G3_results_na.rm$padj<0.001,])
dim(G13_vs_G3_results_na.rm[G13_vs_G3_results_na.rm$padj<0.001,])
dim(G15_vs_G3_results_na.rm[G15_vs_G3_results_na.rm$padj<0.001,])
dim(G17_vs_G3_results_na.rm[G17_vs_G3_results_na.rm$padj<0.001,])
dim(G19_vs_G3_results_na.rm[G19_vs_G3_results_na.rm$padj<0.001,])
dim(G21_vs_G3_results_na.rm[G21_vs_G3_results_na.rm$padj<0.001,])
dim(G23_vs_G3_results_na.rm[G23_vs_G3_results_na.rm$padj<0.001,])
dim(G26_vs_G3_results_na.rm[G26_vs_G3_results_na.rm$padj<0.001,])

# check number of genes DE padj<0.001 and 2x down reg
dim(FLB7_vs_FLB3_results_na.rm[FLB7_vs_FLB3_results_na.rm$padj<0.001 & FLB7_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB10_vs_FLB3_results_na.rm[FLB10_vs_FLB3_results_na.rm$padj<0.001 & FLB10_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB13_vs_FLB3_results_na.rm[FLB13_vs_FLB3_results_na.rm$padj<0.001 & FLB13_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB15_vs_FLB3_results_na.rm[FLB15_vs_FLB3_results_na.rm$padj<0.001 & FLB15_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB17_vs_FLB3_results_na.rm[FLB17_vs_FLB3_results_na.rm$padj<0.001 & FLB17_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB19_vs_FLB3_results_na.rm[FLB19_vs_FLB3_results_na.rm$padj<0.001 & FLB19_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB21_vs_FLB3_results_na.rm[FLB21_vs_FLB3_results_na.rm$padj<0.001 & FLB21_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB23_vs_FLB3_results_na.rm[FLB23_vs_FLB3_results_na.rm$padj<0.001 & FLB23_vs_FLB3_results_na.rm$log2FoldChange<0.5,])
dim(FLB26_vs_FLB3_results_na.rm[FLB26_vs_FLB3_results_na.rm$padj<0.001 & FLB26_vs_FLB3_results_na.rm$log2FoldChange<0.5,])

dim(G7_vs_G3_results_na.rm[G7_vs_G3_results_na.rm$padj<0.001 & G7_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G10_vs_G3_results_na.rm[G10_vs_G3_results_na.rm$padj<0.001 & G10_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G13_vs_G3_results_na.rm[G13_vs_G3_results_na.rm$padj<0.001 & G13_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G15_vs_G3_results_na.rm[G15_vs_G3_results_na.rm$padj<0.001 & G15_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G17_vs_G3_results_na.rm[G17_vs_G3_results_na.rm$padj<0.001 & G17_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G19_vs_G3_results_na.rm[G19_vs_G3_results_na.rm$padj<0.001 & G19_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G21_vs_G3_results_na.rm[G21_vs_G3_results_na.rm$padj<0.001 & G21_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G23_vs_G3_results_na.rm[G23_vs_G3_results_na.rm$padj<0.001 & G23_vs_G3_results_na.rm$log2FoldChange<0.5,])
dim(G26_vs_G3_results_na.rm[G26_vs_G3_results_na.rm$padj<0.001 & G26_vs_G3_results_na.rm$log2FoldChange<0.5,])

# check number of genes DE padj<0.001 and 2x up reg
dim(FLB7_vs_FLB3_results_na.rm[FLB7_vs_FLB3_results_na.rm$padj<0.001 & FLB7_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB10_vs_FLB3_results_na.rm[FLB10_vs_FLB3_results_na.rm$padj<0.001 & FLB10_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB13_vs_FLB3_results_na.rm[FLB13_vs_FLB3_results_na.rm$padj<0.001 & FLB13_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB15_vs_FLB3_results_na.rm[FLB15_vs_FLB3_results_na.rm$padj<0.001 & FLB15_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB17_vs_FLB3_results_na.rm[FLB17_vs_FLB3_results_na.rm$padj<0.001 & FLB17_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB19_vs_FLB3_results_na.rm[FLB19_vs_FLB3_results_na.rm$padj<0.001 & FLB19_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB21_vs_FLB3_results_na.rm[FLB21_vs_FLB3_results_na.rm$padj<0.001 & FLB21_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB23_vs_FLB3_results_na.rm[FLB23_vs_FLB3_results_na.rm$padj<0.001 & FLB23_vs_FLB3_results_na.rm$log2FoldChange>2,])
dim(FLB26_vs_FLB3_results_na.rm[FLB26_vs_FLB3_results_na.rm$padj<0.001 & FLB26_vs_FLB3_results_na.rm$log2FoldChange>2,])

dim(G7_vs_G3_results_na.rm[G7_vs_G3_results_na.rm$padj<0.001 & G7_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G10_vs_G3_results_na.rm[G10_vs_G3_results_na.rm$padj<0.001 & G10_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G13_vs_G3_results_na.rm[G13_vs_G3_results_na.rm$padj<0.001 & G13_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G15_vs_G3_results_na.rm[G15_vs_G3_results_na.rm$padj<0.001 & G15_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G17_vs_G3_results_na.rm[G17_vs_G3_results_na.rm$padj<0.001 & G17_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G19_vs_G3_results_na.rm[G19_vs_G3_results_na.rm$padj<0.001 & G19_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G21_vs_G3_results_na.rm[G21_vs_G3_results_na.rm$padj<0.001 & G21_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G23_vs_G3_results_na.rm[G23_vs_G3_results_na.rm$padj<0.001 & G23_vs_G3_results_na.rm$log2FoldChange>2,])
dim(G26_vs_G3_results_na.rm[G26_vs_G3_results_na.rm$padj<0.001 & G26_vs_G3_results_na.rm$log2FoldChange>2,])




#Are the genes which are differentially expressed in the G at e.g. 21 DAA the same as the ones DE at 23 DAA? 
merged_G23_G21 <- merge(as.data.frame(G23_vs_G3_results_na.rm),as.data.frame(G21_vs_G3_results_na.rm),by="row.names")
head(merged_G23_G21)
dim(merged_G23_G21[merged_G23_G21$padj.x<0.01 & merged_G23_G21$padj.y<0.01,])



bck_res<-bck_res_G10_vs_G3 


# sort results on padj
ordered_res <- bck_res[order(bck_res$padj),]
head(ordered_res)

ordered_res_narm <- na.omit(ordered_res)

ordered_res_filt <- ordered_res_narm[ordered_res_narm$padj<0.05,]

# output ordered_res to csv

write.csv(ordered_res,file="all_vs_FLB3_results.csv")

# make MA plot pdf
pdf(file="MAplot.pdf", height=10, width = 10)
plotMA(bckCDS_1,ylim=c(-2,2),main="DESeq2")
dev.off()

# make MA plot png
png(file="MAplot.png", height=1000, width = 1000)
plotMA(bckCDS_1,ylim=c(-2,2),main="DESeq2")
dev.off()

# make MA plot eps
postscript(file="RNAi_vs_WT_MA.eps", height=100, width = 100)
plotMA(bckCDS_1,ylim=c(-2,2),main="DESeq2")
dev.off()


