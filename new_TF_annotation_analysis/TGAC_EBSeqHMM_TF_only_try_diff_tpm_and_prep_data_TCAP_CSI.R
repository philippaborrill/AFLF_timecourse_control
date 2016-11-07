# Aim is to look to see if differentially expressed TF are expressed over 1 to 10 tpm 
# 23.9.2016


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
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\9_EBSeqHMM_TF_only_new_annotation\\")
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
plot(num_genes_max_tpm, main="Number of diff expr TFs over threshold")


###### TF family info ########
# now want to add in info about TF family
TF_family <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\triticum_aestivum_TFs_final_high_conf.csv", header=TRUE)
head(TF_family)

# add TF family info to tpmData_av
tpmData_av_family <- merge(tpmData_av, TF_family, by.x = "X", by.y= "gene")
head(tpmData_av_family)

##### now do for NAC ####
tpmData_av_NAC <- tpmData_av_family[which(tpmData_av_family$family =="NAC"),]
head(tpmData_av_NAC)
dim(tpmData_av_NAC)

# calculate numbers of genes which are expressed at a max over each value for high conf genes
#initiate empty dataframe
num_genes_high_conf_max_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_high_conf_max_tpm)

for (i in 1:10){
  num_genes <- length(tpmData_av_NAC$maxtpm[which(tpmData_av_NAC$maxtpm >i)])
  num_genes_high_conf_max_tpm <- rbind(num_genes_high_conf_max_tpm, c(i,num_genes))
}
num_genes_high_conf_max_tpm
colnames(num_genes_high_conf_max_tpm) <- num_genes_high_conf_max_tpm[1,]
num_genes_high_conf_max_tpm <- num_genes_high_conf_max_tpm[-1,]
num_genes_high_conf_max_tpm
plot(num_genes_high_conf_max_tpm, main = "Diff expr NACs over tpm threshold")
#######

##### now do for WRKY ####
tpmData_av_WRKY <- tpmData_av_family[which(tpmData_av_family$family =="WRKY"),]
head(tpmData_av_WRKY)
dim(tpmData_av_WRKY)

# calculate numbers of genes which are expressed at a max over each value for high conf genes
#initiate empty dataframe
num_genes_high_conf_max_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_high_conf_max_tpm)

for (i in 1:10){
  num_genes <- length(tpmData_av_WRKY$maxtpm[which(tpmData_av_WRKY$maxtpm >i)])
  num_genes_high_conf_max_tpm <- rbind(num_genes_high_conf_max_tpm, c(i,num_genes))
}
num_genes_high_conf_max_tpm
colnames(num_genes_high_conf_max_tpm) <- num_genes_high_conf_max_tpm[1,]
num_genes_high_conf_max_tpm <- num_genes_high_conf_max_tpm[-1,]
num_genes_high_conf_max_tpm
plot(num_genes_high_conf_max_tpm, main = "Diff expr WRKYs over tpm threshold")
#######

##### now do for MYB ####
tpmData_av_MYB <- tpmData_av_family[which(tpmData_av_family$family =="MYB"),]
head(tpmData_av_MYB)
dim(tpmData_av_MYB)

# calculate numbers of genes which are expressed at a max over each value for high conf genes
#initiate empty dataframe
num_genes_high_conf_max_tpm <- (c("tpm", "num_genes_over_threshold"))
head(num_genes_high_conf_max_tpm)

for (i in 1:10){
  num_genes <- length(tpmData_av_MYB$maxtpm[which(tpmData_av_MYB$maxtpm >i)])
  num_genes_high_conf_max_tpm <- rbind(num_genes_high_conf_max_tpm, c(i,num_genes))
}
num_genes_high_conf_max_tpm
colnames(num_genes_high_conf_max_tpm) <- num_genes_high_conf_max_tpm[1,]
num_genes_high_conf_max_tpm <- num_genes_high_conf_max_tpm[-1,]
num_genes_high_conf_max_tpm
plot(num_genes_high_conf_max_tpm, main = "Diff expr MYBs over tpm threshold")
#######


tpmData_av_conf <- tpmData_av
head(tpmData_av_conf)

#### make file to use with TCAP > 2 tpm ###
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
write.csv(tpmData_2tpm, file="tpmData_av_per_timepoint_diff_expr_TF_2tpm_for_TCAP.csv")

#### make file to use with TCAP > 3 tpm ###
tpmData_3tpm <- tpmData_av_conf[which(tpmData_av_conf$maxtpm>3),]
dim(tpmData_3tpm)
colnames(tpmData_3tpm)
tpmData_3tpm <- tpmData_3tpm[,c(1,4:13)]
head(tpmData_3tpm)
dim(tpmData_3tpm)

# rename colnames to fit with TCAP test data set in cyverse
colnames(tpmData_3tpm) <- c("","3", "7", "10", "13", "15", "17", "19", "21", "23", "26")
head(tpmData_3tpm)
# save as csv
write.csv(tpmData_3tpm, file="tpmData_av_per_timepoint_diff_expr_TF_3tpm_for_TCAP.csv")



# now produce file to use with CSI

# read in counts data (need to log transform to use it)
countData <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\counts_summarised_per_gene.csv")
head(countData)


### do for 2 tpm CSI ###
head(tpmData_2tpm)

# select only genes whcih are differentially expressed with >2 tpm and high confidence (same genes as used for TCAP)
countData_merged <- merge(countData, tpmData_2tpm, by.x="X", by.y="")
head(countData_merged)
dim(countData_merged)

# now only keep FLB count data
countData_merged <- countData_merged[,1:31]
head(countData_merged)

rownames(countData_merged) <- countData_merged$X
countData_merged <- countData_merged[,-1]
head(countData_merged)

head(countData_merged+0.1)

# add 0.1 to each count before log to avoid negative infinity result
logcountData_merge <- log(countData_merged+0.1)
head(logcountData_merge)

# add in correct nomeclature to columns
colnames(logcountData_merge) <- rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3)
head(logcountData_merge)

logcountData_merge2 <- rbind(colnames(logcountData_merge), logcountData_merge)
head(logcountData_merge2)
dim(logcountData_merge2)
colnames(logcountData_merge2) <- rep("WT",30)
head(logcountData_merge2)

write.csv(logcountData_merge2, file="countData_diff_expr_TF_2tpm_high_conf_for_CSI.csv")


# now only select TF from specific TF family (to reduce number of genes in network)
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\9_EBSeqHMM_TF_only_new_annotation\\")
logcountData_merge2 <- read.csv(file="countData_diff_expr_TF_2tpm_high_conf_for_CSI.csv")

# now read in TF family info
TF_info <-  read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\triticum_aestivum_TFs_final_high_conf.csv", header=TRUE)

dim(logcountData_merge2)
dim(TF_info)
head(logcountData_merge2)
head(TF_info)

logcountData_with_TF_family <- merge(logcountData_merge2, TF_info, by.x = "X", by.y ="gene")
head(logcountData_with_TF_family)
dim(logcountData_with_TF_family)


# now select only NAC and WRKY and MYB (to reduce number of genes)
logcountData_NAC_WRKY_MYB <- rbind(logcountData_with_TF_family[which(logcountData_with_TF_family$family =="NAC" ),],
                                   logcountData_with_TF_family[which(logcountData_with_TF_family$family =="WRKY"),],
                                   logcountData_with_TF_family[which(logcountData_with_TF_family$family =="MYB"),])
dim(logcountData_NAC_WRKY_MYB)
head(logcountData_NAC_WRKY_MYB)
tail(logcountData_NAC_WRKY_MYB)

logcountData_NAC_WRKY_MYB <- logcountData_NAC_WRKY_MYB[,1:31]
head(logcountData_NAC_WRKY_MYB)
# add in correct nomeclature to columns
colnames(logcountData_NAC_WRKY_MYB) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(logcountData_NAC_WRKY_MYB)

logcountData_NAC_WRKY_MYB <- rbind(colnames(logcountData_NAC_WRKY_MYB), logcountData_NAC_WRKY_MYB)
head(logcountData_NAC_WRKY_MYB)
dim(logcountData_NAC_WRKY_MYB)
colnames(logcountData_NAC_WRKY_MYB) <- c("condition",rep("WT",30))
head(logcountData_NAC_WRKY_MYB)

write.csv(logcountData_NAC_WRKY_MYB, file="countData_diff_expr_TF_2tpm_high_conf_NAC_WRKY_MYB_for_CSI.csv")

### do for 3 tpm CSI ###
head(tpmData_3tpm)

# select only genes whcih are differentially expressed with >2 tpm and high confidence (same genes as used for TCAP)
countData_merged <- merge(countData, tpmData_3tpm, by.x="X", by.y="")
head(countData_merged)
dim(countData_merged)

# now only keep FLB count data
countData_merged <- countData_merged[,1:31]
head(countData_merged)

rownames(countData_merged) <- countData_merged$X
countData_merged <- countData_merged[,-1]
head(countData_merged)

head(countData_merged+0.1)

# add 0.1 to each count before log to avoid negative infinity result
logcountData_merge <- log(countData_merged+0.1)
head(logcountData_merge)

# add in correct nomeclature to columns
colnames(logcountData_merge) <- rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3)
head(logcountData_merge)

logcountData_merge2 <- rbind(colnames(logcountData_merge), logcountData_merge)
head(logcountData_merge2)
dim(logcountData_merge2)
colnames(logcountData_merge2) <- rep("WT",30)
head(logcountData_merge2)

write.csv(logcountData_merge2, file="countData_diff_expr_TF_3tpm_high_conf_for_CSI.csv")


# now only select TF from specific TF family (to reduce number of genes in network)
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\9_EBSeqHMM_TF_only_new_annotation\\")
logcountData_merge2 <- read.csv(file="countData_diff_expr_TF_3tpm_high_conf_for_CSI.csv")

# now read in TF family info
TF_info <-  read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\triticum_aestivum_TFs_final_high_conf.csv", header=TRUE)

dim(logcountData_merge2)
dim(TF_info)
head(logcountData_merge2)
head(TF_info)

logcountData_with_TF_family <- merge(logcountData_merge2, TF_info, by.x = "X", by.y ="gene")
head(logcountData_with_TF_family)
dim(logcountData_with_TF_family)


# now select only NAC and WRKY and MYB (to reduce number of genes)
logcountData_NAC_WRKY_MYB <- rbind(logcountData_with_TF_family[which(logcountData_with_TF_family$family =="NAC" ),],
                                   logcountData_with_TF_family[which(logcountData_with_TF_family$family =="WRKY"),],
                                   logcountData_with_TF_family[which(logcountData_with_TF_family$family =="MYB"),])
dim(logcountData_NAC_WRKY_MYB)
head(logcountData_NAC_WRKY_MYB)
tail(logcountData_NAC_WRKY_MYB)

logcountData_NAC_WRKY_MYB <- logcountData_NAC_WRKY_MYB[,1:31]
head(logcountData_NAC_WRKY_MYB)
# add in correct nomeclature to columns
colnames(logcountData_NAC_WRKY_MYB) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(logcountData_NAC_WRKY_MYB)

logcountData_NAC_WRKY_MYB <- rbind(colnames(logcountData_NAC_WRKY_MYB), logcountData_NAC_WRKY_MYB)
head(logcountData_NAC_WRKY_MYB)
dim(logcountData_NAC_WRKY_MYB)
colnames(logcountData_NAC_WRKY_MYB) <- c("condition",rep("WT",30))
head(logcountData_NAC_WRKY_MYB)

write.csv(logcountData_NAC_WRKY_MYB, file="countData_diff_expr_TF_3tpm_high_conf_NAC_WRKY_MYB_for_CSI.csv")

##### do for 4 to 10 tpm CSI #####
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\9_EBSeqHMM_TF_only_new_annotation\\")
countData <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\counts_summarised_per_gene.csv")
head(countData)

tpm <- 4

for (tpm in 4:10) {

tpmData_tpm <- tpmData_av_conf[which(tpmData_av_conf$maxtpm>tpm),]
dim(tpmData_tpm)
colnames(tpmData_tpm)
tpmData_tpm <- tpmData_tpm[,c(1,4:13)]
head(tpmData_tpm)



# select only genes whcih are differentially expressed with >(selected level of tpm) and high confidence (same genes as used for TCAP)
countData_merged <- merge(countData, tpmData_tpm, by.x="X", by.y="X")
head(countData_merged)
dim(countData_merged)

# now only keep FLB count data
countData_merged <- countData_merged[,1:31]
head(countData_merged)

rownames(countData_merged) <- countData_merged$X
countData_merged <- countData_merged[,-1]
head(countData_merged)

head(countData_merged+0.1)

# add 0.1 to each count before log to avoid negative infinity result
logcountData_merge <- log(countData_merged+0.1)
head(logcountData_merge)

# add in correct nomeclature to columns
colnames(logcountData_merge) <- rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3)
head(logcountData_merge)

logcountData_merge2 <- rbind(colnames(logcountData_merge), logcountData_merge)
head(logcountData_merge2)
dim(logcountData_merge2)
colnames(logcountData_merge2) <- rep("WT",30)
head(logcountData_merge2)

write.csv(logcountData_merge2, file=paste0("countData_diff_expr_TF_",tpm,"tpm_high_conf_for_CSI.csv"))


# now only select TF from specific TF family (to reduce number of genes in network)
# read in TF family info
TF_info <-  read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\triticum_aestivum_TFs_final_high_conf.csv", header=TRUE)

dim(logcountData_merge2)
dim(TF_info)
head(logcountData_merge2)
head(TF_info)
colnames(logcountData_merge2)

logcountData_with_TF_family <- merge(logcountData_merge2, TF_info, by.x = 0, by.y ="gene")
head(logcountData_with_TF_family)
dim(logcountData_with_TF_family)


# now select only NAC and WRKY and MYB (to reduce number of genes)
logcountData_NAC_WRKY_MYB <- rbind(logcountData_with_TF_family[which(logcountData_with_TF_family$family =="NAC" ),],
                                   logcountData_with_TF_family[which(logcountData_with_TF_family$family =="WRKY"),],
                                   logcountData_with_TF_family[which(logcountData_with_TF_family$family =="MYB"),])
dim(logcountData_NAC_WRKY_MYB)
head(logcountData_NAC_WRKY_MYB)
tail(logcountData_NAC_WRKY_MYB)

logcountData_NAC_WRKY_MYB <- logcountData_NAC_WRKY_MYB[,1:31]
head(logcountData_NAC_WRKY_MYB)
# add in correct nomeclature to columns
colnames(logcountData_NAC_WRKY_MYB) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(logcountData_NAC_WRKY_MYB)

logcountData_NAC_WRKY_MYB <- rbind(colnames(logcountData_NAC_WRKY_MYB), logcountData_NAC_WRKY_MYB)
head(logcountData_NAC_WRKY_MYB)
dim(logcountData_NAC_WRKY_MYB)
colnames(logcountData_NAC_WRKY_MYB) <- c("condition",rep("WT",30))
head(logcountData_NAC_WRKY_MYB)

write.csv(logcountData_NAC_WRKY_MYB, file=paste0("countData_diff_expr_TF_",tpm,"tpm_high_conf_NAC_WRKY_MYB_for_CSI.csv"))

}


# now make file separately for NAC, for WRKY, for MYB at 10 tpm and run to check it doesn't drastically change the network by having 1 family alone
# read in TF family info
TF_info <-  read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\triticum_aestivum_TFs_final_high_conf.csv", header=TRUE)


setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\9_EBSeqHMM_TF_only_new_annotation\\")

count_data_NAC_MYB_WRKY <- read.csv("countData_diff_expr_TF_10tpm_high_conf_NAC_WRKY_MYB_for_CSI.csv", header=T)
head(TF_info)
head(count_data_NAC_MYB_WRKY)
dim(count_data_NAC_MYB_WRKY)

count_data_with_family <- merge(count_data_NAC_MYB_WRKY,TF_info, by.x = "condition", by.y = "gene")
head(count_data_with_family)
dim(count_data_with_family)

# extract just NAC
count_data_NAC <- count_data_with_family[which(count_data_with_family$family == "NAC"),]
dim(count_data_NAC)
head(count_data_NAC)
count_data_NAC <- count_data_NAC[,1:31]

# add in correct nomeclature to columns
colnames(count_data_NAC) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(count_data_NAC)

count_data_NAC <- rbind(colnames(count_data_NAC), count_data_NAC)
head(count_data_NAC)
dim(count_data_NAC)
colnames(count_data_NAC) <- c("condition",rep("WT",30))
head(count_data_NAC)

write.csv(count_data_NAC, file="countData_diff_expr_NAC_10tpm_high_conf_for_CSI.csv", row.names = F)

# extract just MYB
count_data_MYB <- count_data_with_family[which(count_data_with_family$family == "MYB"),]
dim(count_data_MYB)
head(count_data_MYB)
count_data_MYB <- count_data_MYB[,1:31]

# add in correct nomeclature to columns
colnames(count_data_MYB) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(count_data_MYB)

count_data_MYB <- rbind(colnames(count_data_MYB), count_data_MYB)
head(count_data_MYB)
dim(count_data_MYB)
colnames(count_data_MYB) <- c("condition",rep("WT",30))
head(count_data_MYB)

write.csv(count_data_MYB, file="countData_diff_expr_MYB_10tpm_high_conf_for_CSI.csv", row.names = F)

# extract just WRKY
count_data_WRKY <- count_data_with_family[which(count_data_with_family$family == "WRKY"),]
dim(count_data_WRKY)
head(count_data_WRKY)
count_data_WRKY <- count_data_WRKY[,1:31]

# add in correct nomeclature to columns
colnames(count_data_WRKY) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(count_data_WRKY)

count_data_WRKY <- rbind(colnames(count_data_WRKY), count_data_WRKY)
head(count_data_WRKY)
dim(count_data_WRKY)
colnames(count_data_WRKY) <- c("condition",rep("WT",30))
head(count_data_WRKY)

write.csv(count_data_WRKY, file="countData_diff_expr_WRKY_10tpm_high_conf_for_CSI.csv", row.names = F)
