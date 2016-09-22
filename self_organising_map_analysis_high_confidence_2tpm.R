# Aim: to create self organising maps to cluster my gene expression over the time-course
# will use the R package Kohonene which Jemima has been using

### In Jemima's original script she creates a tpm table which has rows as transcripts and each column for a different sample. 
# I want to use the tpm per gene not transcript so I need to use txtimportData to calculate expression per gene rather than per transcript

install.packages("kohonen")

out_dir <- "Y:/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/1_SOM_analysis/2tpm_high_conf/"

#setwd("//nbi//group-data//ifs//NBI//Research-Groups//Cristobal-Uauy//PB_AFLF//control_timecourse//TGAC_kallisto_analysis//kallisto_results_bootstrap//results/")
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results")

source("https://bioconductor.org/biocLite.R")
biocLite("tximportData")
install.packages("readr")
library(tximportData)
library(readr)

# make vector pointing to the kallisto results files   ########
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

# move into directory where I will save this analysis
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis")

# to see counts summarised per gene
head(txi$counts)
colnames(txi$counts)


# save counts summarised per gene
write.csv(txi$counts, file="counts_summarised_per_gene.csv")

# to see tpm summarised per gene
head(txi$abundance)
colnames(txi$abundance)

# save tpm summarised per gene
write.csv(txi$abundance, file="tpm_summarised_per_gene.csv")

# see lengths summarised per gene
head(txi$length)

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="length_per_gene.csv")


#### can restart here to load vector of tpms per gene #####
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis")

tpmData <- read.csv("tpm_summarised_per_gene.csv")
head(tpmData)
colnames(tpmData)
 # adjust rownames
rownames(tpmData) <- tpmData$X
tpmData <- tpmData[,-1]
head(tpmData)

# just use FLB data
tpmData <- tpmData[,1:30]
head(tpmData)
# average per timepoint
tpmData$T3 <- (tpmData[,1] + tpmData[,2] + tpmData[,3]) / 3
tpmData$T7 <- (tpmData[,4] + tpmData[,5] + tpmData[,6]) / 3
tpmData$T10 <- (tpmData[,7] + tpmData[,8] + tpmData[,9]) / 3
tpmData$T13 <- (tpmData[,10] + tpmData[,11] + tpmData[,12]) / 3
tpmData$T15 <- (tpmData[,13] + tpmData[,14] + tpmData[,15]) / 3
tpmData$T17 <- (tpmData[,16] + tpmData[,17] + tpmData[,18]) / 3
tpmData$T19 <- (tpmData[,19] + tpmData[,20] + tpmData[,21]) / 3
tpmData$T21 <- (tpmData[,22] + tpmData[,23] + tpmData[,24]) / 3
tpmData$T23 <- (tpmData[,25] + tpmData[,26] + tpmData[,27]) / 3
tpmData$T26 <- (tpmData[,28] + tpmData[,29] + tpmData[,30]) / 3

head(tpmData)

#create columns saying if expression is less the 2tpm at all time points (want to exclude timepoints less than 2 tpm)
tpmData$test <- ifelse(tpmData$T3 < 2 & tpmData$T7 < 2 & tpmData$T10 < 2 & tpmData$T13 < 2 & tpmData$T15 < 2 & tpmData$T17 < 2 &
                         tpmData$T19 < 2 & tpmData$T21 < 2 & tpmData$T23 < 2 & tpmData$T26 < 2, 0, 1)
head(tpmData$test)
head(tpmData)

# at this stage Jemima also filters to only keep transcripts with GO terms, I won't do this too but might try again later with this filtering criteria
## need to create list of GO terms per gene
######## first generate category mapping (GO information) ##############

genes_GO <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\GO_terms_per_gene.csv")
head(genes_GO)
dim(genes_GO)
colnames(genes_GO)

# need to convert factors to characters to do melt
i <- sapply(genes_GO, is.factor)
genes_GO[i] <- lapply(genes_GO[i], as.character)
head(genes_GO)
dim(genes_GO)
str(genes_GO)

#melt to have 1 gene per line with 1 GO term
library(reshape2)
genes_GO2 <- melt(genes_GO, id.vars=c("Gene_ID"))
head(genes_GO2,100)
tail(genes_GO2,100)
genes_GO2 <- melt(genes_GO, id=c("genes_GO$Gene_ID" ))
head(genes_GO2)

colnames(genes_GO2) <- c("Gene","GO_group","GO_term")
head(genes_GO2,20)
dim(genes_GO2)

# select only rows with a GO_term (it kept GO1 to GO61 for all transcripts even though most don't have GO10 +)

melted_gene_info_complete <- genes_GO2[ grep("GO",genes_GO2$GO_term) ,]
head(melted_gene_info_complete)
dim(melted_gene_info_complete)

# save csv file
write.csv(melted_gene_info_complete, file="Y:\\PB_AFLF\\control_timecourse\\GO_terms_per_gene_long_form.csv")

##### read in csv file of GO per gene (can use again)#####
melted_gene_info_complete2 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\GO_terms_per_gene_long_form.csv")
head(melted_gene_info_complete2)
melted_gene_info_complete2 <- melted_gene_info_complete2[,2:4]

# make dataframe to use in GOseq
GO_data.frame <- as.data.frame(cbind(as.character(melted_gene_info_complete2$Gene),melted_gene_info_complete2$GO_term))
colnames(GO_data.frame) <- c("id", "GO")
head(GO_data.frame)
dim(GO_data.frame)
str(GO_data.frame)
is.data.frame(GO_data.frame)

#since there are many GO terms for one gene, lots of rows get duplicated, so for now remove GO terms and duplicates and add them back later
head(tpmData)
tpmData$target_id <- rownames(tpmData)
head(tpmData)
colnames(tpmData)
dim(tpmData)
tpmF <- tpmData[, c(42, 1:41)]
rownames(tpmF) <- NULL
tpmF <- unique(tpmF)
dim(tpmF)
head(tpmF)

############ignoring GO filtering #############
head(tpmData)
tpmData$target_id <- rownames(tpmData)
head(tpmData)
colnames(tpmData)
tpmF <- tpmData[, c(42, 1:41)]
rownames(tpmF) <- NULL
dim(tpmF)
head(tpmF)

#filter out genes that have average tpm less than 2 at all time points
tpmF2 <- tpmF[tpmF$test == 1, ]
dim(tpmF2)
head(tpmF2)

# now only keep high confidence genes
# now want to do the same for high confidence genes therefore need to get list of high confidence genes
gene_conf <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_confidence.csv", header=TRUE)
head(gene_conf)
colnames(gene_conf) <- c("gene","confidence")
head(gene_conf)
# add confidence level to tpmF2
tpmF2 <- merge(tpmF2, gene_conf, by.x = "target_id", by.y= "gene")
tpmF2 <- tpmF2[tpmF2$confidence == "High", ]
dim(tpmF2)
head(tpmF2)


######now filter for only DE genes between subsequent time-points ########### 
# later might want to try genes which are differentially expressed between 3 DAA and all other timepoints
# or doing it for all genes (i.e not just diff expr)
#load the tables that have the DEseq output

FLB7_FLB3 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB7 _vs_ FLB3 _results.csv", header=TRUE)
FLB10_FLB7 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB10 _vs_ FLB7 _results.csv", header=TRUE)
FLB13_FLB10 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB13 _vs_ FLB10 _results.csv", header=TRUE)
FLB15_FLB13 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB15 _vs_ FLB13 _results.csv", header=TRUE)
FLB17_FLB15 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB17 _vs_ FLB15 _results.csv", header=TRUE)
FLB19_FLB17 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB19 _vs_ FLB17 _results.csv", header=TRUE)
FLB21_FLB19 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB21 _vs_ FLB19 _results.csv", header=TRUE)
FLB23_FLB21 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB23 _vs_ FLB21 _results.csv", header=TRUE)
FLB26_FLB23 <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\DESeq2_analysis\\FLB26 _vs_ FLB23 _results.csv", header=TRUE)

head(FLB26_FLB23)

# filter for qvalue 0.01 and then rbind all DE gene lists together to get final list of DE genes
FLB7_FLB3 <- FLB7_FLB3[FLB7_FLB3$padj < 0.01,]
FLB10_FLB7 <- FLB10_FLB7[FLB10_FLB7$padj < 0.01,]
FLB13_FLB10 <- FLB13_FLB10[FLB13_FLB10$padj < 0.01,]
FLB15_FLB13 <- FLB15_FLB13[FLB15_FLB13$padj < 0.01,]
FLB17_FLB15 <- FLB17_FLB15[FLB17_FLB15$padj < 0.01,]
FLB19_FLB17 <- FLB19_FLB17[FLB19_FLB17$padj < 0.01,]
FLB21_FLB19 <- FLB21_FLB19[FLB21_FLB19$padj < 0.01,]
FLB23_FLB21 <- FLB23_FLB21[FLB23_FLB21$padj < 0.01,]
FLB26_FLB23 <- FLB26_FLB23[FLB26_FLB23$padj < 0.01,]

DE_table <- rbind(FLB7_FLB3, FLB10_FLB7, FLB13_FLB10, FLB15_FLB13, FLB17_FLB15, FLB19_FLB17, FLB21_FLB19, FLB23_FLB21, FLB26_FLB23)
DE_genes <- unique(DE_table$X)
head(DE_table)
head(DE_genes)
length(DE_genes)

##### check if NAM genes are in DE_genes list ##########
grep("TRIAE_CS42_6AS_TGACv1_486738_AA1564640",DE_genes)
grep("TRIAE_CS42_6BS_TGACv1_513229_AA1635270",DE_genes)
grep("TRIAE_CS42_6DS_TGACv1_542626_AA1725630",DE_genes)
grep("TRIAE_CS42_2AS_TGACv1_113243_AA0353410",DE_genes)
grep("TRIAE_CS42_2BS_TGACv1_145996_AA0452240",DE_table)
grep("TRIAE_CS42_2DS_TGACv1_179582_AA0608070",DE_genes)

# only 2AS copy is in the list ....


#### if only want to use DE genes with 2tpm #######
head(tpmF2)
tpmF3 <- subset(tpmF2, target_id %in% DE_genes)
write.csv(tpmF3, file = "tpmData_DEonly.csv", row.names = TRUE, quote = FALSE)
head(tpmF3)
dim(tpmF3)

#now proceed with just average tpm across reps for now
colnames(tpmF3)
colnames(tpmF3)[c(1,32:41)]
av_tpm <- tpmF3[, c(1,32:41)]
rownames(av_tpm) <- av_tpm[,1]
av_tpm[,1] <- NULL
av_tpm <- as.matrix(av_tpm)
head(av_tpm)

##### if want to use 2 tpm but don't care if DE between adjacent timepoints #########
#now proceed with just average tpm across reps for now
colnames(tpmF2)
colnames(tpmF2)[c(1,32:41)]
av_tpm <- tpmF2[, c(1,32:41)]
rownames(av_tpm) <- av_tpm[,1]
av_tpm[,1] <- NULL
av_tpm <- as.matrix(av_tpm)
head(av_tpm)
dim(av_tpm)

#scale by centering the mean at 0 and the sd to 1 per gene
means <- apply(av_tpm, 1, mean)
sds <- apply(av_tpm, 1, sd)
av_tpm_sc <- (av_tpm - means) / sds
head(av_tpm_sc)

# ASK JEMIMA what this does
test <- complete.cases(av_tpm_sc[, c(1:10)])
not <- which(test == FALSE)
if (unique(test) == TRUE) {
  print("No FALSE")
} else {
  av_tpm_sc <- av_tpm_sc[-not,]
}


#### now test out kohonen #######
library("kohonen")
#set a seed (randomly chosen number) so that the results are reproducible by run as it is a random process
set.seed(42)


## ask Jemima to explain this part
#use eigen values of the input matrix to define an aspect ratio for the self organising map
total_number_of_nodes <- 100
ratio <- sort(svd(as.matrix(av_tpm_sc))$d, decreasing=TRUE)[1] / sort(svd(as.matrix(av_tpm_sc))$d, decreasing=TRUE)[2]
grid_y <- sqrt(total_number_of_nodes / ratio)
grid_y <- c(floor(grid_y), ceiling(grid_y))
grid_y <- grid_y[which(grid_y%%2==0)]
grid_x <- floor(ratio * grid_y)
grid_y
grid_x
#this gives x = 12, y = 6 for DE genes padj <0.01, tpm > 2
# for all genes expressed over 2 tpm gives x =13, y=8
# for all genes > 2 tpm high conf gives x =14, y =8

# set working dir
setwd("Y:/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/1_SOM_analysis/2tpm_high_conf")


# plot the map of line graphs of expression trends
av_tpm_sc.som <- som(data = av_tpm_sc, grid = somgrid(grid_x, grid_y, "hexagonal"), toroidal = TRUE)
map <- plot(av_tpm_sc.som, main = "FLB", codeRendering = "lines")
tiff(filename = "FLB_SOM_all_genes.tiff",width = 1000, height = 800, units = "px")
plot(av_tpm_sc.som, main = "FLB", codeRendering = "lines")
dev.off()

#check that enough iterations are being done
graphics.off()
plot(av_tpm_sc.som, main = "FLB", type = "changes")
tiff(filename = "FLB_SOM_changes_all_genes.tiff")
plot(av_tpm_sc.som, main = "FLB", type = "changes")
dev.off()

# check that you have similar numbers of genes in each node (https://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/)
plot(av_tpm_sc.som, type="count")

# check distance between neighbours
plot(av_tpm_sc.som, type="dist.neighbours")

# number of groups
groups <- 13

tiff(filename = paste0(out_dir, "FLB_SOM_hclust_all_genes_",groups,"_groups.tiff", sep = ""), width = 1200, height = 500)
plot(hclust(dist(av_tpm_sc.som$codes)), cex=0.8)
rect.hclust(hclust(dist(av_tpm_sc.som$codes)),k=groups, border = "blue")
dev.off()


#you can add boundaries showing which nodes would cluster together if you clustered them based on hierachical clustering, but again you have to decide on the number of groups eg 10
som.hc <- cutree(hclust(dist(av_tpm_sc.som$codes)), groups)
plot(av_tpm_sc.som, main = "FLB", codeRendering = "lines")
add.cluster.boundaries(av_tpm_sc.som, som.hc)

tiff(filename = paste0(out_dir, "FLB_SOM_",groups,"_groups_all_transcript.tiff", sep = ""), width = 1000, height = 500)
som.hc <- cutree(hclust(dist(av_tpm_sc.som$codes)), groups)
plot(av_tpm_sc.som, main = "FLB", codeRendering = "lines")
add.cluster.boundaries(av_tpm_sc.som, som.hc)
dev.off()


library(RColorBrewer)
# try adding background colour according to group - r color brewer with legend
tiff(filename = paste0(out_dir, "FLB_SOM_",groups,"_groups_all_transcript_coloured_brewer_palette_outline_legend.tiff", sep = ""), width = 1000, height = 500)
som.hc <- cutree(hclust(dist(av_tpm_sc.som$codes)), groups)
plot(av_tpm_sc.som, main = "FLB", codeRendering = "lines",  palette.name=palette(brewer.pal(12,"Set3")), bgcol=som.hc, ncolors=groups )
legend(-1.05,6, legend=unique(som.hc), fill=palette(brewer.pal(12,"Set3")), bty="n")
add.cluster.boundaries(av_tpm_sc.som, som.hc, lwd=3)
dev.off()

# try adding background colour according to group - r color brewer with legend
tiff(filename = paste0(out_dir, "FLB_SOM_",groups,"_groups_all_transcript_coloured_brewer_palette_no_outline_legend.tiff", sep = ""), width = 1000, height = 500)
som.hc <- cutree(hclust(dist(av_tpm_sc.som$codes)), groups)
plot(av_tpm_sc.som, main = "FLB", codeRendering = "lines",  palette.name=palette(brewer.pal(12,"Set3")), bgcol=som.hc, ncolors=groups )
legend(-1.05,6, legend=unique(som.hc), fill=palette(brewer.pal(12,"Set3")), bty="n")
dev.off()


#you can get access to the different elements of the kohonen object using $...the unit.classif is a vector with cluster 
names(av_tpm_sc.som)
clusters <- av_tpm_sc.som$data
clusters <- cbind(clusters, av_tpm_sc.som$unit.classif)
head(clusters)
colnames(clusters)
colnames(clusters) <- c(colnames(clusters)[1:10], "cluster")
head(clusters)
groups <- cbind(seq(1, length(som.hc)), som.hc)
colnames(groups) <- c("cluster", "hier_group")
groups <- data.frame(groups)
groups <- cbind(groups, av_tpm_sc.som$codes)
head(groups)

write.table(clusters, file = paste0(out_dir, "FLB_clusters_11.tsv", sep = ""), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(groups, file = paste0(out_dir, "FLB_codes_groups_11.tsv", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#### read in information about lengths and GO terms for each cluster #########
# start inspecting each cluster (this will save each hierachical group of genes and then perform GO term enrichment analysis using GOseq on each group
clusters <- data.frame(clusters)

# read in GO terms
all_go <- read.csv("Y:\\PB_AFLF\\control_timecourse\\GO_terms_per_gene_long_form.csv")
head(all_go)
all_go <- all_go[,c(2,4)]
head(all_go)


#create vector for gene_lengths

# need to get lengths of genes not of transcripts
lengths <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\length_per_gene.csv", header=T)
head(lengths)
colnames(lengths) <- c("gene", "length")
head(lengths)

t1 <- subset(lengths, gene %in% rownames(av_tpm_sc))
head(t1)
dim(t1)
# turn into a vector called gene.lens to use with GOSeq
gene.lens <- as.numeric(t1$length)
names(gene.lens) = t1$gene
head(gene.lens)


####### Do GO term enrichment and plot graph for each group ####

assayed.genes <- as.vector(rownames(av_tpm_sc))
length(assayed.genes)

library(goseq)
#i=1

for (i in seq(1, length(unique(groups$hier_group)))) {
	group_i <- groups[groups$hier_group == i, 1]
	groupi_genes <- subset(clusters, cluster %in% group_i)
	groupi_genes$target_id <- rownames(groupi_genes)
	write.table(groupi_genes, file = paste0(out_dir, "group_", i, "_genes.tsv", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE)
	#get the GO terms
	## Philippa thinks this step is unnecessary groupi_GO <- merge(groupi_genes, all_go, all.x = TRUE, all.y = FALSE)
	#now do GO stats analysis on the cluster using the cluster genes as the defined set and all genes used for clustering as the gene universe??
	#create a named binary vector for genes where one means expressed and 0 means not expressed
	de.genes <- groupi_genes$target_id
	gene.vector=as.integer(assayed.genes%in%de.genes)
	names(gene.vector)=assayed.genes
	head(gene.vector)
	#now carry out the GOseq analysis
	pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
	GO.wall = goseq(pwf, gene2cat = all_go)
	#this gave table with p-values...now correct for multiple testing using FDR
	enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]
	head(enriched.GO)
	# add new column with over represented GO terms padj
	GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")
	# add new column with under represented GO terms padj
	GO.wall$under_rep_padj=p.adjust(GO.wall$under_represented_pvalue, method="BH")
	dep.GO=GO.wall$category[p.adjust(GO.wall$under_represented_pvalue, method="BH")<.05]
	write.table(GO.wall, file = paste0(out_dir, "group_", i, "_GOseq.tsv", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE)
	
}

#now plot the behaviour of each of the clusters
library("reshape2")
library("ggplot2")

i=1

for (i in seq(1, length(unique(groups$hier_group)))) {
	file_path <- paste0(out_dir, "group_", i, "_genes.tsv", sep = "")
	groupi_genes <- read.table(file_path, sep = "\t", header = TRUE)
	groupi <- groupi_genes[, c(12, 1:10)]
	rownames(groupi) <- NULL
	colnames(groupi) <- c("Time_point", colnames(clusters)[1:10])
	groupi <- data.frame(t(groupi))
	groupi <- groupi[c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11), ]
	groupi$Time_point <- c("T03","T07","T10","T13","T15","T17","T19","T21","T23","T26")
	x <- ncol(groupi)
	y <- x-1
	groupi <- groupi[, c(x, 1:y)]
	groupi.melted <- melt(groupi, id = "Time_point")
	groupi.melted$value <- as.numeric(groupi.melted$value)
	
	#define the values for the mean of the group (all transcripts in group averaged)
	T03 <- mean(groupi_genes$T3)
	T07 <- mean(groupi_genes$T7)
	T10 <- mean(groupi_genes$T10)
	T13 <- mean(groupi_genes$T13)
	T15 <- mean(groupi_genes$T15)	
	T17 <- mean(groupi_genes$T17)
	T19 <- mean(groupi_genes$T19)
	T21 <- mean(groupi_genes$T21)
	T23 <- mean(groupi_genes$T23)
	T26 <- mean(groupi_genes$T26)
	groupi.means <- data.frame(T03,T07,T10,T13,T15,T17,T19,T21,T23,T26)
	head(groupi.means)
	groupi.means <- data.frame(t(groupi.means))
	groupi.means$Time_point <- rownames(groupi.means)
	groupi.means <- groupi.means[, c(2, 1)]
	colnames(groupi.means) <- c("Time_point", "mean")
	means.melted <- melt(groupi.means, id = "Time_point")
	means.melted$value <- as.numeric(means.melted$value)
	
	# define the mean of the clusters within the group (e.g. group 1 has 5 clusters included, what is the mean of them. 
	# NB. each cluster is already the mean of the genes inside that cluster)
	hier_i <- groups[groups$hier_group == i,]
	T03 <- mean(hier_i$T3)
	T07 <- mean(hier_i$T7)
	T10 <- mean(hier_i$T10)
	T13 <- mean(hier_i$T13)
	T15 <- mean(hier_i$T15)	
	T17 <- mean(hier_i$T17)
	T19 <- mean(hier_i$T19)
	T21 <- mean(hier_i$T21)
	T23 <- mean(hier_i$T23)
	T26 <- mean(hier_i$T26)
	hier.means <- data.frame(T03,T07,T10,T13,T15,T17,T19,T21,T23,T26)
	hier.means <- data.frame(t(hier.means))
	hier.means$Time_point <- rownames(hier.means)
	hier.means <- hier.means[, c(2, 1)]
	colnames(hier.means) <- c("Time_point", "mean")
	hmeans.melted <- melt(hier.means, id = "Time_point")
	hmeans.melted$value <- as.numeric(hmeans.melted$value)
	
	n_transcripts <- nrow(groupi_genes)
	
	plot_title <- paste0("FLB_group ", i, ": all transcripts")
	plot_name <- paste0(out_dir, "FLB_group_", i, "_genes.tiff", sep = "")
	
	#initiate the plot
	plot <- ggplot(groupi.melted, aes(x = Time_point, y = value))
	#add the underlay lines
	plot <- plot + geom_line(data = groupi.melted, mapping = aes(x = Time_point, y = value, group = variable), colour = "LIGHTGRAY", alpha = 1/2, size = 1/2)
	#add the mean line
	plot <- plot + geom_line(data = means.melted, mapping = aes(x = Time_point, y = value, group = variable, color = "Gene mean"), size=1)
	plot <- plot + geom_line(data = hmeans.melted, mapping = aes(x = Time_point, y = value, group = variable, color = "Cluster mean") , size=1)
	plot <- plot + scale_color_manual(values = c("Cluster mean" = "RED", "Gene mean" = "BLACK"))
	plot <- plot + ggtitle(plot_title) + labs( y = "Normalised tpm") + annotate("text", label = paste0("n = ", n_transcripts), x= Inf , hjust = 1, y=Inf, vjust = 1)
	ggsave(plot, file = plot_name)
}


##### want to try plotting the behaviour of the cluster with the NAM genes in it, with the NAM gene expression shown too #####
# NAM genes are in cluster 5

out_dir <- "Y:/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/1_SOM_analysis/2tpm/2nd_run/"

i= 5
file_path <- paste0(out_dir, "group_", i, "_genes.tsv", sep = "")
groupi_genes <- read.table(file_path, sep = "\t", header = TRUE)
groupi <- groupi_genes[, c(12, 1:10)]
rownames(groupi) <- NULL
colnames(groupi) <- c("Time_point", colnames(clusters)[1:10])
groupi <- data.frame(t(groupi))
groupi <- groupi[c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11), ]
groupi$Time_point <- c("T03","T07","T10","T13","T15","T17","T19","T21","T23","T26")
x <- ncol(groupi)
y <- x-1
groupi <- groupi[, c(x, 1:y)]
groupi.melted <- melt(groupi, id = "Time_point")
groupi.melted$value <- as.numeric(groupi.melted$value)

#define the values for the mean of the group (all transcripts in group averaged)
T03 <- mean(groupi_genes$T3)
T07 <- mean(groupi_genes$T7)
T10 <- mean(groupi_genes$T10)
T13 <- mean(groupi_genes$T13)
T15 <- mean(groupi_genes$T15)	
T17 <- mean(groupi_genes$T17)
T19 <- mean(groupi_genes$T19)
T21 <- mean(groupi_genes$T21)
T23 <- mean(groupi_genes$T23)
T26 <- mean(groupi_genes$T26)
groupi.means <- data.frame(T03,T07,T10,T13,T15,T17,T19,T21,T23,T26)
head(groupi.means)
groupi.means <- data.frame(t(groupi.means))
groupi.means$Time_point <- rownames(groupi.means)
groupi.means <- groupi.means[, c(2, 1)]
colnames(groupi.means) <- c("Time_point", "mean")
means.melted <- melt(groupi.means, id = "Time_point")
means.melted$value <- as.numeric(means.melted$value)

# define the mean of the clusters within the group (e.g. group 1 has 5 clusters included, what is the mean of them. 
# NB. each cluster is already the mean of the genes inside that cluster)
hier_i <- groups[groups$hier_group == i,]
T03 <- mean(hier_i$T3)
T07 <- mean(hier_i$T7)
T10 <- mean(hier_i$T10)
T13 <- mean(hier_i$T13)
T15 <- mean(hier_i$T15)	
T17 <- mean(hier_i$T17)
T19 <- mean(hier_i$T19)
T21 <- mean(hier_i$T21)
T23 <- mean(hier_i$T23)
T26 <- mean(hier_i$T26)
hier.means <- data.frame(T03,T07,T10,T13,T15,T17,T19,T21,T23,T26)
hier.means <- data.frame(t(hier.means))
hier.means$Time_point <- rownames(hier.means)
hier.means <- hier.means[, c(2, 1)]
colnames(hier.means) <- c("Time_point", "mean")
hmeans.melted <- melt(hier.means, id = "Time_point")
hmeans.melted$value <- as.numeric(hmeans.melted$value)


# now get the NAM gene expression

file_path <- paste0(out_dir, "group_", i, "_genes.tsv", sep = "")
groupi_genes <- read.table(file_path, sep = "\t", header = TRUE)
groupi <- groupi_genes[, c(12, 1:10)]
head(groupi)
groupi_tr <- data.frame(t(groupi))
head(groupi_tr)

groupi_tr <- data.frame(groupi_tr$TRIAE_CS42_6AS_TGACv1_486738_AA1564640, groupi_tr$TRIAE_CS42_6DS_TGACv1_542626_AA1725630,
                   groupi_tr$TRIAE_CS42_2AS_TGACv1_113243_AA0353410, groupi_tr$TRIAE_CS42_2BS_TGACv1_145996_AA0452240, 
                   groupi_tr$TRIAE_CS42_2DS_TGACv1_179582_AA0608070)
head(groupi_tr)
colnames(groupi_tr)
colnames(groupi_tr) <- gsub("groupi_tr.","",colnames(groupi_tr))
colnames(groupi_tr)
groupi_tr <- groupi_tr[2:11,]
head(groupi_tr)
# add timepoint column
groupi_tr$Time_point <- hmeans.melted$Time_point
head(groupi_tr)
rownames(groupi_tr) <- NULL
head(groupi_tr)

# convert factors to numeric
indx <- sapply(groupi_tr, is.factor)
groupi_tr[indx] <- lapply(groupi_tr[indx], function(x) as.numeric(as.character(x)))

is.factor(groupi_tr$TRIAE_CS42_6AS_TGACv1_486738_AA1564640)
is.numeric(groupi_tr$TRIAE_CS42_6AS_TGACv1_486738_AA1564640)
head(groupi_tr)
melted_groupi_tr <- melt(groupi_tr)
head(melted_groupi_tr)

n_transcripts <- nrow(groupi_genes)

plot_title <- paste0("FLB_group ", i, ": all transcripts_with_NAM")
plot_name <- paste0(out_dir, "FLB_group_", i, "_genes_with_NAM.tiff", sep = "")

#initiate the plot
plot <- ggplot(groupi.melted, aes(x = Time_point, y = value))
#add the underlay lines
plot <- plot + geom_line(data = groupi.melted, mapping = aes(x = Time_point, y = value, group = variable), colour = "LIGHTGRAY", alpha = 1/2, size = 1/2)
#add the mean line
plot <- plot + geom_line(data = means.melted, mapping = aes(x = Time_point, y = value, group = variable, color = "Gene mean"), size=1)
plot <- plot + geom_line(data = hmeans.melted, mapping = aes(x = Time_point, y = value, group = variable, color = "Cluster mean") , size=1)
plot <- plot + geom_line(data = melted_groupi_tr, mapping = aes(x = Time_point, y = value, group=variable, color = variable), size=1)
plot <- plot + scale_color_manual(values = c("Cluster mean" = "RED", "Gene mean" = "BLACK", "TRIAE_CS42_6AS_TGACv1_486738_AA1564640" = "yellow", "TRIAE_CS42_6DS_TGACv1_542626_AA1725630"="chartreuse3", "TRIAE_CS42_2AS_TGACv1_113243_AA0353410" = "cornflowerblue", "TRIAE_CS42_2BS_TGACv1_145996_AA0452240"= "blue", "TRIAE_CS42_2DS_TGACv1_179582_AA0608070" = "dark green"))
plot <- plot + ggtitle(plot_title) + labs( y = "Normalised tpm") + annotate("text", label = paste0("n = ", n_transcripts), x= Inf , hjust = 1, y=Inf, vjust = 1)
ggsave(plot, file = plot_name)
