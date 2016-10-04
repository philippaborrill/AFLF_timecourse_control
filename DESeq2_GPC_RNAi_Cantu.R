#Run DESeq2
#Philippa Borrill
#27-11-13 edited 3/6/15 edited 03/10/2016


# or set WD for using PC version
setwd("Y:\\PB_AFLF\\control_timecourse\\GPC_RNAi_Cantu\\per_gene")

library("DESeq2")
install.packages("DeSeq2")

bckCountTable <- read.table("counts_per_gene_GPC_RNAi_Cantu.txt", header=TRUE, row.names=1)
head(bckCountTable)
colnames(bckCountTable)
colnames(bckCountTable) <- c("RNAi_rep1", "RNAi_rep2", "RNAi_rep3","RNAi_rep4","WT_rep1","WT_rep2", "WT_rep3")
head(bckCountTable)
row.names(bckCountTable)

# get rid of any rows which have a sum of 0
counts <- bckCountTable[rowSums(bckCountTable) != 0,]
dim(counts)
head(counts)
dim(bckCountTable)
samples <- data.frame(row.names=c("RNAi_rep1","RNAi_rep2","RNAi_rep3","RNAi_rep4","WT_rep1","WT_rep2", "WT_rep3"), condition=as.factor(c(rep("RNAi",4),rep("WT",3))))
samples

# round counts to integers (required for DESeq)
counts <- round(counts)
head(counts)

bckCDS <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
bckCDS

#make sure the WT is used and the reference condition:
colData(bckCDS)$condition <- relevel(colData(bckCDS)$condition, "WT")

bckCDS_1 <- DESeq(bckCDS)

bck_res <- results(bckCDS_1)

head(bck_res)

write.csv(bck_res,file="RNAi_vs_WT_results.csv")


pdf(file="RNAi_vs_WT_MA.pdf", height=10, width = 10)
plotMA(bckCDS_1,ylim=c(-2,2),main="DESeq2")
dev.off()

mcols(bck_res, use.names=TRUE)



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
install.packages("vsn")
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

pdf(file="deseq2_heatmaps_samplebysample.pdf", height =8, width =8)
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



