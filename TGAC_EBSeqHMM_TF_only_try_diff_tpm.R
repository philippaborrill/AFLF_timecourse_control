# Aim is to look to see if differentially expressed genes are expressed over 2 and 3 tpm 



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

### now read in differentially expressed gene list from EBSeqHMM ####
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\4b_EBSeq_HMM_TF_only\\")
DE_genes <- read.csv(file="GeneDECalls_FLB_FDR0.001.csv", header=T)
head(DE_genes)
dim(DE_genes)

# merge together DE_genes with tpmData_max
DE_genes <- merge(DE_genes,tpmData_av, by.x = "X", by.y = 0)
head(DE_genes)

tpmData_av <- DE_genes
head(tpmData_av)
dim(tpmData_av)

#####all genes######
# calculate numbers of genes which are expressed at a max over each value
#initiate empty dataframe
num_genes_max_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_max_tpm)

for (i in 1:10){
  num_genes <- length(tpmData_av$maxtpm[which(tpmData_av$maxtpm >i)])
  num_genes_max_tpm <- rbind(num_genes_max_tpm, c(i,num_genes))
}
num_genes_max_tpm
colnames(num_genes_max_tpm) <- num_genes_max_tpm[1,]
num_genes_max_tpm <- num_genes_max_tpm[-1,]
num_genes_max_tpm
plot.new()
par(mfrow=c(1,1))
plot(num_genes_max_tpm)

###### high confidence genes ########
# now want to do the same for high confidence genes therefore need to get list of high confidence genes
gene_conf <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_confidence.csv", header=TRUE)
head(gene_conf)
colnames(gene_conf) <- c("gene","confidence")
head(gene_conf)
# add confidence level to tpmData_av
tpmData_av_conf <- merge(tpmData_av, gene_conf, by.x = "X", by.y= "gene")
head(tpmData_av_conf)

dim(tpmData_av)
dim(gene_conf)
dim(tpmData_av_conf)

tpmData_av_conf <- tpmData_av_conf[which(tpmData_av_conf$confidence =="High"),]
head(tpmData_av_conf)
dim(tpmData_av_conf)

# calculate numbers of genes which are expressed at a max over each value for high conf genes
#initiate empty dataframe
num_genes_high_conf_max_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_high_conf_max_tpm)

for (i in 1:10){
  num_genes <- length(tpmData_av_conf$maxtpm[which(tpmData_av_conf$maxtpm >i)])
  num_genes_high_conf_max_tpm <- rbind(num_genes_high_conf_max_tpm, c(i,num_genes))
}
num_genes_high_conf_max_tpm
colnames(num_genes_high_conf_max_tpm) <- num_genes_high_conf_max_tpm[1,]
num_genes_high_conf_max_tpm <- num_genes_high_conf_max_tpm[-1,]
num_genes_high_conf_max_tpm
plot(num_genes_high_conf_max_tpm)

# save file of high confidence genes with > 2 tpm for use with TCAP

tpmData_high_conf_2tpm <- tpmData_av_conf[which(tpmData_av_conf$maxtpm>2),]
dim(tpmData_high_conf_2tpm)
colnames(tpmData_high_conf_2tpm)
tpmData_high_conf_2tpm <- tpmData_high_conf_2tpm[,c(1,4:13)]
head(tpmData_high_conf_2tpm)
dim(tpmData_high_conf_2tpm)

# rename colnames to fit with TCAP test data set in cyverse
colnames(tpmData_high_conf_2tpm) <- c("","3", "7", "10", "13", "15", "17", "19", "21", "23", "26")
head(tpmData_high_conf_2tpm)

# save as csv
write.csv(tpmData_high_conf_2tpm, file="tpmData_av_per_timepoint_diff_expr_TF_2tpm_high_conf.csv")

# Want to use all genes (not just high confidence)
tpmData_av_conf <- merge(tpmData_av, gene_conf, by.x = "X", by.y= "gene")
head(tpmData_av_conf)

tpmData_2tpm <- tpmData_av_conf[which(tpmData_av_conf$maxtpm>2),]
dim(tpmData_2tpm)
colnames(tpmData_2tpm)
tpmData_2tpm <- tpmData_2tpm[,c(1,4:13)]
head(tpmData_2tpm)
dim(tpmData_2tpm)

# rename colnames to fit with TCAP test data set in cyverse
colnames(tpmData_2tpm) <- c("","3", "7", "10", "13", "15", "17", "19", "21", "23", "26")
head(tpmData_2tpm)
# save as csv
write.csv(tpmData_2tpm, file="tpmData_av_per_timepoint_diff_expr_TF_2tpm.csv")



# now produce file to use with CSI

# read in counts data (need to log transform to use it)
countData <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\counts_summarised_per_gene.csv")
head(countData)
head(tpmData_high_conf_2tpm)

# select only genes whcih are differentially expressed with >2 tpm and high confidence (same genes as used for TCAP)
countData_merged <- merge(countData, tpmData_high_conf_2tpm, by.x="X", by.y="")
head(countData_merged)

# now only keep FLB count data
countData_merged <- countData_merged[,1:31]
head(countData_merged)

rownames(countData_merged) <- countData_merged$X
countData_merged <- countData_merged[,-1]
head(countData_merged)

logcountData_merge <- log(countData_merged)
head(logcountData_merge)

# add in correct nomeclature to columns
colnames(logcountData_merge) <- rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3)
head(logcountData_merge)

logcountData_merge2 <- rbind(colnames(logcountData_merge), logcountData_merge)
head(logcountData_merge2)
dim(logcountData_merge2)
colnames(logcountData_merge2) <- rep("WT",30)
head(logcountData_merge2)

#write.csv(logcountData_merge2, file="countData_diff_expr_TF_2tpm_high_conf.csv")
write.csv(logcountData_merge2, file="countData_diff_expr_TF_2tpm_high_conf_repeat.csv")

# now only select TF from specific TF family (to reduce number of genes in network)
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\4b_EBSeq_HMM_TF_only\\")
logcountData_merge2 <- read.csv(file="countData_diff_expr_TF_2tpm_high_conf.csv")

# now read in TF family info
TF_info <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_TF_family.csv")

dim(logcountData_merge2)
dim(TF_info)
head(logcountData_merge2)
head(TF_info)

logcountData_with_TF_family <- merge(logcountData_merge2, TF_info, by.x = "condition", by.y ="X.2.Gene_ID")
head(logcountData_with_TF_family)
dim(logcountData_with_TF_family)


# now select only NAC and WRKY (to reduce number of genes)
logcountData_NAC_WRKY <- rbind(logcountData_with_TF_family[which(logcountData_with_TF_family$X.19.TF_family =="NAC" ),],
                               logcountData_with_TF_family[which(logcountData_with_TF_family$X.19.TF_family =="WRKY" ),])
dim(logcountData_NAC_WRKY)
head(logcountData_NAC_WRKY)
tail(logcountData_NAC_WRKY)

logcountData_NAC_WRKY <- logcountData_NAC_WRKY[,1:31]
head(logcountData_NAC_WRKY)
# add in correct nomeclature to columns
colnames(logcountData_NAC_WRKY) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(logcountData_NAC_WRKY)

logcountData_NAC_WRKY <- rbind(colnames(logcountData_NAC_WRKY), logcountData_NAC_WRKY)
head(logcountData_NAC_WRKY)
dim(logcountData_NAC_WRKY)
colnames(logcountData_NAC_WRKY) <- c("condition",rep("WT",30))
head(logcountData_NAC_WRKY)

write.csv(logcountData_NAC_WRKY, file="countData_diff_expr_TF_2tpm_high_conf_NAC_WRKY.csv")

