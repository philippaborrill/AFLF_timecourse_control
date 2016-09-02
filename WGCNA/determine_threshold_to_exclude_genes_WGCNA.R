# Determine tpm threshold to exclude genes from wGCNA #
# 02.09.2016 #
# Philippa Borrill # 

#### can restart here to load vector of tpms per gene #####
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\")
options(scipen=10000)
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

# just keep the average per timepoint
head(tpmData)
colnames(tpmData)[31:40]
tpmData_av <- tpmData[,31:40]
head(tpmData_av)

tpmData_av$maxtpm <- apply(tpmData_av[,1:10],1,max)
head(tpmData_av)

tpmData_av$avtpm <- apply(tpmData_av[,1:10],1,mean)
head(tpmData_av)

#####all genes######
# calculate numbers of genes which are expressed at a max over each value
#initiate empty dataframe
num_genes_max_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_max_tpm)

for (i in -1:10){
num_genes <- length(tpmData_av$maxtpm[which(tpmData_av$maxtpm >i)])
num_genes_max_tpm <- rbind(num_genes_max_tpm, c(i,num_genes))
}
num_genes_max_tpm
colnames(num_genes_max_tpm) <- num_genes_max_tpm[1,]
num_genes_max_tpm <- num_genes_max_tpm[-1,]
num_genes_max_tpm
plot(num_genes_max_tpm)

###### high confidence genes ########
# now want to do the same for high confidence genes therefore need to get list of high confidence genes
gene_conf <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_confidence.csv", header=TRUE)
head(gene_conf)
colnames(gene_conf) <- c("gene","confidence")
head(gene_conf)
# add confidence level to tpmData_av
tpmData_av_conf <- merge(tpmData_av, gene_conf, by.x = 0, by.y= "gene")
head(tpmData_av_conf)

dim(tpmData_av)
dim(gene_conf)
dim(tpmData_av_conf)

tpmData_av_conf <- tpmData_av_conf[which(tpmData_av_conf$confidence =="High"),]
head(tpmData_av_conf)

# calculate numbers of genes which are expressed at a max over each value for high conf genes
#initiate empty dataframe
num_genes_high_conf_max_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_high_conf_max_tpm)

for (i in -1:10){
  num_genes <- length(tpmData_av_conf$maxtpm[which(tpmData_av_conf$maxtpm >i)])
  num_genes_high_conf_max_tpm <- rbind(num_genes_high_conf_max_tpm, c(i,num_genes))
}
num_genes_high_conf_max_tpm
colnames(num_genes_high_conf_max_tpm) <- num_genes_high_conf_max_tpm[1,]
num_genes_high_conf_max_tpm <- num_genes_high_conf_max_tpm[-1,]
num_genes_high_conf_max_tpm
plot(num_genes_high_conf_max_tpm)

###### transcription factors ########
# now want to do the same for transcription factor genes
# load file which only contains genes which are TFs
TF_family <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_TF_family.csv", header =T)
head(TF_family)
colnames(TF_family) <- c("gene", "TF_family")
head(TF_family)

# select only genes from tpmData_av which are TFs
tpmData_av_TF <- merge(tpmData_av, TF_family, by.x = 0, by.y= "gene")
head(tpmData_av_TF)
dim(tpmData_av_TF)
dim(TF_family)

# calculate numbers of genes which are expressed at a max over each value for high conf genes
#initiate empty dataframe
num_genes_TF_max_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_TF_max_tpm)

for (i in -1:10){
  num_genes <- length(tpmData_av_TF$maxtpm[which(tpmData_av_TF$maxtpm >i)])
  num_genes_TF_max_tpm <- rbind(num_genes_TF_max_tpm, c(i,num_genes))
}
num_genes_TF_max_tpm
colnames(num_genes_TF_max_tpm) <- num_genes_TF_max_tpm[1,]
num_genes_TF_max_tpm <- num_genes_TF_max_tpm[-1,]
num_genes_TF_max_tpm
plot(num_genes_TF_max_tpm)

##### now combine the 3 sets (all genes, high confidence, and TF) and plot ########

merged_max_data <- merge(num_genes_max_tpm, num_genes_high_conf_max_tpm, by="tpm")
merged_max_data

merged_max_data <- merge(merged_max_data, num_genes_TF_max_tpm, by="tpm")
merged_max_data
colnames(merged_max_data) <- c("tpm", "all_genes", "high_confidence_genes", "transcription_factors")
merged_max_data <- merged_max_data[order(merged_max_data$all_genes, decreasing=T),]
merged_max_data

library(reshape2)
melted_merged_max_data <- melt(merged_max_data, id.vars =c("tpm"))
melted_merged_max_data

summary(melted_merged_max_data)
melted_merged_max_data$tpm <- as.numeric(melted_merged_max_data$tpm)
melted_merged_max_data$value <- as.numeric(melted_merged_max_data$value)

library(ggplot2)
p <- ggplot(data=melted_merged_max_data, aes(x=tpm, y=value, colour=variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("purple", "hotpink", "#666666", "darkblue")) +
    scale_x_continuous(breaks = -1:10) + 
  ylab("Number of genes") +
  scale_fill_manual(guide = FALSE)


p

# put each graph on separate axis to see TF better
p2 <- p + facet_grid( variable ~ ., scales="free_y")

out_dir <- "Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\2_Leaf_WGCNA\\"

library(gridExtra)

pdf(file=paste0(out_dir,"Number_of_genes_with_max_expression_exceeding_tpm_levels.pdf"), height=8, width =9)
grid.arrange(p, p2, nrow=2, ncol=1, top=textGrob("Number of genes with max gene expression exceeding tpm -1 to 10", gp=gpar(fontsize=15)))
dev.off()



####### for average tpm #######

#####all genes######
# calculate numbers of genes which are expressed at a max over each value
#initiate empty dataframe
num_genes_av_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_av_tpm)

for (i in -1:10){
  num_genes <- length(tpmData_av$avtpm[which(tpmData_av$avtpm >i)])
  num_genes_av_tpm <- rbind(num_genes_av_tpm, c(i,num_genes))
}
num_genes_av_tpm
colnames(num_genes_av_tpm) <- num_genes_av_tpm[1,]
num_genes_av_tpm <- num_genes_av_tpm[-1,]
num_genes_av_tpm
plot(num_genes_av_tpm)

###### high confidence genes ########
# now want to do the same for high confidence genes therefore need to get list of high confidence genes
gene_conf <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_confidence.csv", header=TRUE)
head(gene_conf)
colnames(gene_conf) <- c("gene","confidence")
head(gene_conf)
# add confidence level to tpmData_av
tpmData_av_conf <- merge(tpmData_av, gene_conf, by.x = 0, by.y= "gene")
head(tpmData_av_conf)

dim(tpmData_av)
dim(gene_conf)
dim(tpmData_av_conf)

tpmData_av_conf <- tpmData_av_conf[which(tpmData_av_conf$confidence =="High"),]
head(tpmData_av_conf)

# calculate numbers of genes which are expressed at a av over each value for high conf genes
#initiate empty dataframe
num_genes_high_conf_av_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_high_conf_av_tpm)

for (i in -1:10){
  num_genes <- length(tpmData_av_conf$avtpm[which(tpmData_av_conf$avtpm >i)])
  num_genes_high_conf_av_tpm <- rbind(num_genes_high_conf_av_tpm, c(i,num_genes))
}
num_genes_high_conf_av_tpm
colnames(num_genes_high_conf_av_tpm) <- num_genes_high_conf_av_tpm[1,]
num_genes_high_conf_av_tpm <- num_genes_high_conf_av_tpm[-1,]
num_genes_high_conf_av_tpm
plot(num_genes_high_conf_av_tpm)

###### transcription factors ########
# now want to do the same for transcription factor genes
# load file which only contains genes which are TFs
TF_family <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_TF_family.csv", header =T)
head(TF_family)
colnames(TF_family) <- c("gene", "TF_family")
head(TF_family)

# select only genes from tpmData_av which are TFs
tpmData_av_TF <- merge(tpmData_av, TF_family, by.x = 0, by.y= "gene")
head(tpmData_av_TF)
dim(tpmData_av_TF)
dim(TF_family)

# calculate numbers of genes which are expressed at a av over each value for high conf genes
#initiate empty dataframe
num_genes_TF_av_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_TF_av_tpm)

for (i in -1:10){
  num_genes <- length(tpmData_av_TF$avtpm[which(tpmData_av_TF$avtpm >i)])
  num_genes_TF_av_tpm <- rbind(num_genes_TF_av_tpm, c(i,num_genes))
}
num_genes_TF_av_tpm
colnames(num_genes_TF_av_tpm) <- num_genes_TF_av_tpm[1,]
num_genes_TF_av_tpm <- num_genes_TF_av_tpm[-1,]
num_genes_TF_av_tpm
plot(num_genes_TF_av_tpm)

##### now combine the 3 sets (all genes, high confidence, and TF) and plot ########

merged_av_data <- merge(num_genes_av_tpm, num_genes_high_conf_av_tpm, by="tpm")
merged_av_data

merged_av_data <- merge(merged_av_data, num_genes_TF_av_tpm, by="tpm")
merged_av_data
colnames(merged_av_data) <- c("tpm", "all_genes", "high_confidence_genes", "transcription_factors")
merged_av_data <- merged_av_data[order(merged_av_data$all_genes, decreasing=T),]
merged_av_data

library(reshape2)
melted_merged_av_data <- melt(merged_av_data, id.vars =c("tpm"))
melted_merged_av_data

summary(melted_merged_av_data)
melted_merged_av_data$tpm <- as.numeric(melted_merged_av_data$tpm)
melted_merged_av_data$value <- as.numeric(melted_merged_av_data$value)

library(ggplot2)
p <- ggplot(data=melted_merged_av_data, aes(x=tpm, y=value, colour=variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("purple", "hotpink", "#666666", "darkblue")) +
  scale_x_continuous(breaks = -1:10) + 
  ylab("Number of genes") +
  scale_fill_manual(guide = FALSE)


p

# put each graph on separate axis to see TF better
p2 <- p + facet_grid( variable ~ ., scales="free_y")

out_dir <- "Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\2_Leaf_WGCNA\\"

library(gridExtra)

pdf(file=paste0(out_dir,"Number_of_genes_with_mean_expression_exceeding_tpm_levels.pdf"), height=8, width =9)
grid.arrange(p, p2, nrow=2, ncol=1, top=textGrob("Number of genes with mean gene expression exceeding tpm -1 to 10", gp=gpar(fontsize=15)))
dev.off()

# get NAM gene expression

tpmData_av

