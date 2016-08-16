#Run DESeq2
#Philippa Borrill
#27-11-13 edited 3/6/15 edited 09/02/2016
# 


#source R-3.0.2
#R


library("DESeq2")

setwd("/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/")

#setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results")


count.data<-read.table(file="final_output_counts.txt",header=T)  

head(count.data)
dim(count.data)

# make rownames correct
rownames(count.data) <- count.data[,1]
counts <- count.data[,-1]
head(counts)
head(row.names(counts))

# round counts to integers (required for DESeq)
counts <- round(counts)
# get rid of any rows which have a sum of 0
counts <- counts[rowSums(counts) != 0,]
# check this has removed some rows
head(row.names(counts))

# how many genes were expressed:
dim(counts)
head(colnames(counts))

# how many genes expressed in FLB or G separately?
counts_FLB <- counts[,1:30]
head(counts_FLB)
counts_FLB <- counts_FLB[rowSums(counts_FLB) != 0,]
dim(counts_FLB)

counts_G <- counts[,31:60]
head(counts_G)
counts_G <- counts_G[rowSums(counts_G) != 0,]
dim(counts_G)

# assign levels to columns to use as factors

timepoints <- c("FLB3", "FLB7", "FLB10", "FLB13", "FLB15", "FLB17", "FLB19", "FLB21", "FLB23","FLB26","G3", "G7", "G10", "G13", "G15", "G17", "G19", "G21", "G23","G26")
CondVector <- rep(timepoints,each=3)
CondVector

samples <- data.frame(row.names=colnames(counts), condition=as.factor(CondVector))
samples


# get data into DESeq form
bckCDS <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
bckCDS

#make sure the FLB3 is used as the reference condition :
colData(bckCDS)$condition <- relevel(colData(bckCDS)$condition, "FLB3")


# run DeSeq2
bckCDS_1 <- DESeq(bckCDS)

dds <- bckCDS_1

#try out using http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/ tutorial:

#transform the raw counts so we can do clustering:
rld <- rlogTransformation(dds, blind=TRUE)

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

#plot graph to see effect of transformation
pdf (file="DESeq2_VST_and_log2.pdf", height=5, width = 5)
par(mai=ifelse(1:5 <= 2, par("mai"), 0))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright", legend = c(expression("variance stabilizing transformation"), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

#plot graph to show comparison of transformation methods
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

pdf(file="DESeq2_heatmap2reglog.pdf", height=8, width = 8)
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

pdf(file="DESeq2_heatmap3varstab.pdf", height=8, width = 8)
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

#calc sample to sample distances to make dendrogram to look at clustering

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))

pdf(file="deseq2_heatmaps_samplebysample_rld.pdf", height =8, width =8)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

# think I can separate out the FLB and G from the matrix
head(mat)
dim(mat)
FLB_mat<- mat[1:30,1:30]
dim(FLB_mat)

pdf(file="deseq2_heatmaps_samplebysample_rld_FLB.pdf", height =8, width =8)
heatmap.2(FLB_mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

G_mat<- mat[31:60,31:60]
dim(G_mat)

pdf(file="deseq2_heatmaps_samplebysample_rld_G.pdf", height =8, width =8)
heatmap.2(G_mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()


#try PCA 
pdf(file="deseq2_pca_rld.pdf", height =8, width =10)
print(plotPCA(vsd, intgroup=c("condition")))

dev.off()

# 1st separate rld
FLB_rld <- rld[,1:30]
pdf(file="deseq2_pca_rld_FLB.pdf", height =8, width =10)
print(plotPCA(FLB_rld, intgroup=c("condition")))

dev.off()

# 1st separate rld
G_rld <- rld[,31:60]
pdf(file="deseq2_pca_rld_G.pdf", height =8, width =10)
print(plotPCA(G_rld, intgroup=c("condition")))

dev.off()





#try with vsd (rather than rld) to look at clustering # NB this actually doesn't look nice- the WT and RNAi no longer cluster separately! ignore?!

distsVSD <- dist(t(assay(vsd)))
mat <- as.matrix(distsVSD)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))



pdf(file="deseq2_heatmaps_samplebysample_VSD.pdf", height =8, width =8)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
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



# get data for timepoint vs. previous timepoint
i <- "G26"
j <- "G23"

bck_res <- results(bckCDS_1,contrast=c("condition",i,j))

# sort results on padj
ordered_res <- bck_res[order(bck_res$padj),]
head(ordered_res)

# output ordered_res to csv
write.csv(ordered_res,file=paste(i,"_vs_",j,"_results.csv"))
ordered_res_na_rm <- na.omit(ordered_res)

# number of genes padj < 0.05
dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.05,])

# number of genes padj < 0.01
dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01,])

# check number of genes DE padj<0.01 and 2x down reg
dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01 & ordered_res_na_rm$log2FoldChange<0.5,])

# check number of genes DE padj<0.01 and 2x up reg
dim(ordered_res_na_rm[ordered_res_na_rm$padj<0.01 & ordered_res_na_rm$log2FoldChange>2,])



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
pdf(file="RNAi_vs_WT_MA.pdf", height=10, width = 10)
plotMA(bckCDS_1,ylim=c(-2,2),main="DESeq2")
dev.off()

# make MA plot png
png(file="RNAi_vs_WT_MA.png", height=1000, width = 1000)
plotMA(bckCDS_1,ylim=c(-2,2),main="DESeq2")
dev.off()

# make MA plot eps
postscript(file="RNAi_vs_WT_MA.eps", height=100, width = 100)
plotMA(bckCDS_1,ylim=c(-2,2),main="DESeq2")
dev.off()

#output metadata to csv
write.csv((mcols(bck_res, use.names=TRUE)),file="RNAi_vs_WT_metadata.csv")

#assign bck_CDS1 to variable name dds to fit with manual

dds <- bckCDS_1

#try out using http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/ tutorial:

#transform the raw counts so we can do clustering:
rld <- rlogTransformation(dds, blind=TRUE)

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

#plot graph to see effect of transformation
pdf (file="DESeq2_VST_and_log2.pdf", height=5, width = 5)
par(mai=ifelse(1:5 <= 2, par("mai"), 0))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright", legend = c(expression("variance stabilizing transformation"), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

#plot graph to show comparison of transformation methods
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

pdf(file="DESeq2_heatmap2reglog.pdf", height=8, width = 8)
heatmap.2(assay(rld)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale="none",
dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

pdf(file="DESeq2_heatmap3varstab.pdf", height=8, width = 8)
heatmap.2(assay(vsd)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale="none",
dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

#calc sample to sample distances to make dendrogram to look at clustering

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
paste(condition, type, sep=" : "))

pdf(file="deseq2_heatmaps_samplebysample_rld.pdf", height =8, width =8)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

#try with vsd (rather than rld) to look at clustering # NB this actually doesn't look nice- the WT and RNAi no longer cluster separately! ignore?!

distsVSD <- dist(t(assay(vsd)))
mat <- as.matrix(distsVSD)
rownames(mat) <- colnames(mat) <- with(colData(dds),
paste(condition, type, sep=" : "))

pdf(file="deseq2_heatmaps_samplebysample_VSD.pdf", height =8, width =8)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

#try PCA 
pdf(file="deseq2_pca.pdf", height =8, width =10)
print(plotPCA(rld, intgroup=c("condition")))

dev.off()

#trying out replacing outliers with Trimmed mean approach substitutions


ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < .1,
cleaned = results(ddsClean)$padj < .1)
addmargins(tab)
write.csv(as.data.frame(tab),file="sim_condition_treated_results_cleaned_summary_deseq2.csv")
resClean <- results(ddsClean)
write.csv(as.data.frame(resClean),file="sim_condition_treated_results_cleaned_deseq2.csv")


# plot dispersions
pdf (file="deseq2_plot_dispersions.pdf", height = 8, width = 8)
plotDispEsts(dds)
dev.off()

#remove any tests which have little chance of pass to avoid multiple testing error. filtering threshold
attr(bck_res,"filterThreshold")

pdf(file="deseq2_filtering_threshold.pdf", height =8, width=8)

plot(attr(bck_res,"filterNumRej"),type="b", ylab="number of rejections")

dev.off()


# Plot of the maximum Cook’s distance per gene over the rank of the Wald statistics for the condition.  This isn't a result it's just a filtering criteria (already implemented and useful to look at to check the filtering was appropriate.

pdf(file="deseq2_cooksdist.pdf", height=8, width=8)
W <- mcols(dds)$WaldStatistic_condition_RNAi_vs_WT
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)

plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",
ylab="maximum Cook's distance per gene",
ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))

dev.off()

#more independent filtering- plot -log10 of pvalue against mean of normalized counts

pdf(file="deseq2_indep_filt.pdf", height = 8, width = 8)
plot(bck_res$baseMean+1, -log10(bck_res$pvalue),
log="x", xlab="mean of normalized counts",
ylab=expression(-log[10](pvalue)),
ylim=c(0,30),
cex=.4, col=rgb(0,0,0,.3))

dev.off()

#graph of very low count values which will be removed to reduce chance of making a type I error:

pdf (file= "deseq2_remove_low_counts_filt.pdf", height=8, width =8)
use <- bck_res$baseMean > attr(bck_res,"filterThreshold")
table(use)
h1 <- hist(bck_res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(bck_res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c('do not pass'="khaki", 'pass'="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))
dev.off()

#graph of very low count values which will be removed to reduce chance of making a type I error part 2:
# This graph on ranks the p-values from smallest to biggest (x-axis) and plots them. 
# The black line is the actual p-value numbers. 
# The red line has a slope that represents the number of tests divided by the false discovery rate (0.1). 
# The idea here is the FDR is controlled at the 0.1% value for all tests that occur to the left of the right-most intersection of the black and red line.

resFilt <- bck_res[use & !is.na(bck_res$pvalue),]
orderInPlot <- order(resFilt$pvalue)
showInPlot <- (resFilt$pvalue[orderInPlot] <= 0.08)
alpha <- 0.1

pdf (file= "deseq2_remove_low_counts_filt2.pdf", height=8, width =8)
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)
dev.off()

################################ignore below this
#making graphs doesn't work

library(gplots)
library(RColorBrewer)
library(Biobase)
library(oligo)

SizeFactors_bckCDS <- estimateSizeFactors (bckCDS)

Disp_bckCDS <- estimateDispersions(SizeFactors_bckCDS)

vsd_bck <-  varianceStabilizingTransformation(SizeFactors_bckCDS)

pdf(file="DeSeq_heatmap_top30.pdf", height=10, width = 7)
select = order(rowMeans(counts(SizeFactors_bckCDS)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsd_bck)[select,], col = hmcol, trace="none", margin=c(10, 6))
dev.off()


#############



#plot graphs of NAM expression

tpm.data<-read.table(file="final_output_tpm.txt",header=T)  

head(tpm.data)
dim(tpm.data)

# make rownames correct
rownames(tpm.data) <- tpm.data[,1]
tpm <- tpm.data[,-1]
head(tpm)
head(row.names(tpm))

## plot graph sophie genes 
gene9 <- "Traes_5AS_E5FD82E36"
gene10 <- "Traes_5BS_9F043AE17"

### plot graphs
jpeg(filename="Gene9 and gene10 expression tpm flag leaf.jpg", width=480, height = 480, unit="px")
plot.new()
par(mfrow=c(5,1))
par(mar=c(1,4,1,1))
par(oma=c(7,0,0,0))
par(cex.lab=1.5, cex.axis =1)
colours = c(rep("red",3), rep("orange",3), rep("yellow",3), rep("lightgreen",3), rep("green",3), rep("cyan",3), rep("turquoise",3),  rep("blue",3), rep("violet",3), rep("pink",3) )

gene9<- tpm["Traes_5AS_E5FD82E36",1:30]
barplot(as.matrix(gene9), las=2, xaxt="n", ylim=c(0,30), ylab="tpm", col=colours)
legend("topright",legend="gene9",box.lty=0,cex=1)

gene10 <- tpm["Traes_5BS_9F043AE17",1:30]
barplot(as.matrix(gene10), las=2,  ylim=c(0,30), ylab="tpm")
legend("topright",legend="gene10",box.lty=0,cex=1)

dev.off()



### plot graphs
jpeg(filename="NAM expression tpm all tissues.jpg", width=480, height = 480, unit="px")
plot.new()
par(mfrow=c(5,1))
par(mar=c(1,4,1,1))
par(oma=c(7,0,0,0))
par(cex.lab=1.5, cex.axis =1)
colours = c(rep("red",3), rep("orange",3), rep("yellow",3), rep("lightgreen",3), rep("green",3), rep("cyan",3), rep("turquoise",3),  rep("blue",3), rep("violet",3), rep("pink",3) )

NAMA1<- tpm["Traes_6AS_6F89CC969",]
barplot(as.matrix(NAMA1), las=2, xaxt="n", ylim=c(0,100), ylab="tpm", col=colours)
legend("topright",legend="NAM-A1",box.lty=0,cex=2)

NAMD1 <- tpm["Traes_6DS_90A627C3F",]
barplot(as.matrix(NAMD1), las=2, xaxt="n", ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-D1",box.lty=0,cex=2)

NAMA2 <- tpm["Traes_2AS_09457EDD8",]
barplot(as.matrix(NAMA2), las=2, xaxt="n", ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-A2",box.lty=0,cex=2)

NAMB2 <- tpm["Traes_2BS_58AA829F4",]
barplot(as.matrix(NAMB2), las=2, xaxt="n", ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-B2",box.lty=0,cex=2)

NAMD2 <- tpm["Traes_2DS_5DF921ABB",]
barplot(as.matrix(NAMD2), las=2, ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-D2",box.lty=0,cex=2)

dev.off()

## just FLB
jpeg(filename="NAM expression tpm FLB.jpg", width=480, height = 480, unit="px")
plot.new()
par(mfrow=c(5,1))
par(mar=c(1,4,1,1))
par(oma=c(7,0,0,0))
par(cex.lab=1.5, cex.axis =1.2)

NAMA1<- tpm["Traes_6AS_6F89CC969",1:30]
barplot(as.matrix(NAMA1), las=2, xaxt="n", ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-A1",box.lty=0,cex=2)

NAMD1 <- tpm["Traes_6DS_90A627C3F",1:30]
barplot(as.matrix(NAMD1), las=2, xaxt="n", ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-D1",box.lty=0,cex=2)

NAMA2 <- tpm["Traes_2AS_09457EDD8",1:30]
barplot(as.matrix(NAMA2), las=2, xaxt="n", ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-A2",box.lty=0,cex=2)

NAMB2 <- tpm["Traes_2BS_58AA829F4",1:30]
barplot(as.matrix(NAMB2), las=2, xaxt="n", ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-B2",box.lty=0,cex=2)

NAMD2 <- tpm["Traes_2DS_5DF921ABB",1:30]
barplot(as.matrix(NAMD2), las=2, ylim=c(0,100), ylab="tpm")
legend("topright",legend="NAM-D2",box.lty=0,cex=2)

dev.off()

## just G
jpeg(filename="NAM expression tpm G.jpg", width=480, height = 480, unit="px")
plot.new()
par(mfrow=c(5,1))
par(mar=c(1,4,1,1))
par(oma=c(7,0,0,0))
par(cex.lab=1.5, cex.axis =1.2)

NAMA1<- tpm["Traes_6AS_6F89CC969",31:60]
barplot(as.matrix(NAMA1), las=2, xaxt="n", ylim=c(0,10), ylab="tpm")
legend("topright",legend="NAM-A1",box.lty=0,cex=2)

NAMD1 <- tpm["Traes_6DS_90A627C3F",31:60]
barplot(as.matrix(NAMD1), las=2, xaxt="n", ylim=c(0,10), ylab="tpm")
legend("topright",legend="NAM-D1",box.lty=0,cex=2)

NAMA2 <- tpm["Traes_2AS_09457EDD8",31:60]
barplot(as.matrix(NAMA2), las=2, xaxt="n", ylim=c(0,10), ylab="tpm")
legend("topright",legend="NAM-A2",box.lty=0,cex=2)

NAMB2 <- tpm["Traes_2BS_58AA829F4",31:60]
barplot(as.matrix(NAMB2), las=2, xaxt="n", ylim=c(0,10), ylab="tpm")
legend("topright",legend="NAM-B2",box.lty=0,cex=2)

NAMD2 <- tpm["Traes_2DS_5DF921ABB",31:60]
barplot(as.matrix(NAMD2), las=2, ylim=c(0,10), ylab="tpm")
legend("topright",legend="NAM-D2",box.lty=0,cex=2)

dev.off()
