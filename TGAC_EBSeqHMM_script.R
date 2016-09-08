# Script to run EBSeqHMM on cluster to look at genes differentially expressed across time in RNA-seq timecourse
# Philippa Borrill 08.02.2016


# load packages
library("EBSeq")
library("EBSeqHMM")


# using the Vignette as a guide: http://bioconductor.org/packages/release/bioc/vignettes/EBSeqHMM/inst/doc/EBSeqHMM_vignette.pdf
# requires matrix of counts, with samples as columns, individual genes as rows
# the object "conditions" should be a factor with the same length as the number of samples e.g. 1 level per sample
# e.g. if you have 5 timepoints with 3 reps you would have t1,t1,t1,t2,t2,t2 etc
# I think I will need to analys the FLB and G samples separately

## read in data and only select genes with > 1 tpm ##

setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/1_SOM_analysis/")

count.data<-read.csv(file="counts_summarised_per_gene.csv",header=T)  

head(count.data)
dim(count.data)

# make rownames correct
rownames(count.data) <- count.data[,1]
counts <- count.data[,-1]
head(counts)

# now select only genes with > 1tpm in at least 1 timepoint FLB

tpmData_FLB <- read.csv(file="tpm_summarised_per_gene.csv", header=T)
head(tpmData_FLB)
colnames(tpmData_FLB)
# adjust rownames
rownames(tpmData_FLB) <- tpmData_FLB$X
tpmData_FLB <- tpmData_FLB[,-1]
head(tpmData_FLB)

# just use FLB data
tpmData_FLB <- tpmData_FLB[,1:30]
head(tpmData_FLB)


# average per timepoint
tpmData_FLB$T3 <- (tpmData_FLB[,1] + tpmData_FLB[,2] + tpmData_FLB[,3]) / 3
tpmData_FLB$T7 <- (tpmData_FLB[,4] + tpmData_FLB[,5] + tpmData_FLB[,6]) / 3
tpmData_FLB$T10 <- (tpmData_FLB[,7] + tpmData_FLB[,8] + tpmData_FLB[,9]) / 3
tpmData_FLB$T13 <- (tpmData_FLB[,10] + tpmData_FLB[,11] + tpmData_FLB[,12]) / 3
tpmData_FLB$T15 <- (tpmData_FLB[,13] + tpmData_FLB[,14] + tpmData_FLB[,15]) / 3
tpmData_FLB$T17 <- (tpmData_FLB[,16] + tpmData_FLB[,17] + tpmData_FLB[,18]) / 3
tpmData_FLB$T19 <- (tpmData_FLB[,19] + tpmData_FLB[,20] + tpmData_FLB[,21]) / 3
tpmData_FLB$T21 <- (tpmData_FLB[,22] + tpmData_FLB[,23] + tpmData_FLB[,24]) / 3
tpmData_FLB$T23 <- (tpmData_FLB[,25] + tpmData_FLB[,26] + tpmData_FLB[,27]) / 3
tpmData_FLB$T26 <- (tpmData_FLB[,28] + tpmData_FLB[,29] + tpmData_FLB[,30]) / 3

# just keep the average per timepoint
head(tpmData_FLB)
colnames(tpmData_FLB)[31:40]
tpmData_FLB_av <- tpmData_FLB[,31:40]
head(tpmData_FLB_av)

tpmData_FLB_av$maxtpm <- apply(tpmData_FLB_av[,1:10],1,max)
head(tpmData_FLB_av)

# clean up workspace to remove unnecessary dataframes
rm(tpmData_FLB)

# merge together counts_conf  and tpmData_FLB_av
counts_max_FLB <- merge(counts, tpmData_FLB_av, by.x = 0, by.y = 0)
head(counts_max_FLB)
dim(counts_max_FLB)
# select only rows with a maxtpm >1
counts_FLB <- counts_max_FLB[which(counts_max_FLB$maxtpm>1),]

# make rownames correct
rownames(counts_FLB) <- counts_FLB[,1]
counts_FLB <- counts_FLB[,-1]
head(counts_FLB)
head(row.names(counts_FLB))
# remove tpm and maxtpm columns
counts_FLB <- counts_FLB[,1:30]
head(counts_FLB)
dim(counts_FLB)

#convert to matrix
count_matrix_FLB <- as.matrix(counts_FLB)
print("dimensions of count_matrix all genes FLB >1 tpm:")
dim(count_matrix_FLB)

# rename to fit rest of script
FLB_count_matrix <- count_matrix_FLB


# now select only genes with > 1tpm in at least 1 timepoint G

tpmData_G <- read.csv(file="tpm_summarised_per_gene.csv", header=T)
head(tpmData_G)
colnames(tpmData_G)
# adjust rownames
rownames(tpmData_G) <- tpmData_G$X
tpmData_G <- tpmData_G[,-1]
head(tpmData_G)

# just use G data
tpmData_G <- tpmData_G[,31:60]
head(tpmData_G)


# average per timepoint
tpmData_G$T3 <- (tpmData_G[,1] + tpmData_G[,2] + tpmData_G[,3]) / 3
tpmData_G$T7 <- (tpmData_G[,4] + tpmData_G[,5] + tpmData_G[,6]) / 3
tpmData_G$T10 <- (tpmData_G[,7] + tpmData_G[,8] + tpmData_G[,9]) / 3
tpmData_G$T13 <- (tpmData_G[,10] + tpmData_G[,11] + tpmData_G[,12]) / 3
tpmData_G$T15 <- (tpmData_G[,13] + tpmData_G[,14] + tpmData_G[,15]) / 3
tpmData_G$T17 <- (tpmData_G[,16] + tpmData_G[,17] + tpmData_G[,18]) / 3
tpmData_G$T19 <- (tpmData_G[,19] + tpmData_G[,20] + tpmData_G[,21]) / 3
tpmData_G$T21 <- (tpmData_G[,22] + tpmData_G[,23] + tpmData_G[,24]) / 3
tpmData_G$T23 <- (tpmData_G[,25] + tpmData_G[,26] + tpmData_G[,27]) / 3
tpmData_G$T26 <- (tpmData_G[,28] + tpmData_G[,29] + tpmData_G[,30]) / 3

# just keep the average per timepoint
head(tpmData_G)
colnames(tpmData_G)[31:40]
tpmData_G_av <- tpmData_G[,31:40]
head(tpmData_G_av)

tpmData_G_av$maxtpm <- apply(tpmData_G_av[,1:10],1,max)
head(tpmData_G_av)

# clean up workspace to remove unnecessary dataframes
rm(tpmData_G)

# merge together counts  and tpmData_G_av
counts_max_G <- merge(counts, tpmData_G_av, by.x = 0, by.y = 0)
head(counts_max_G)
dim(counts_max_G)
# select only rows with a maxtpm >1
counts_G <- counts_max_G[which(counts_max_G$maxtpm>1),]

# make rownames correct
rownames(counts_G) <- counts_G[,1]
counts_G <- counts_G[,-1]
head(counts_G)
head(row.names(counts_G))
# remove tpm and maxtpm columns
counts_G <- counts_G[,1:30]
head(counts_G)
dim(counts_G)

#convert to matrix
count_matrix_G <- as.matrix(counts_G)
print("dimensions of count_matrix all genes G >1 tpm:")
dim(count_matrix_G)

# rename to fit rest of script
G_count_matrix <- count_matrix_G


### set up vectors of conditions #####
# make vector of time-points
timepoints <- c("t3", "t7", "t10", "t13", "t15", "t17", "t19", "t21", "t23","t26")
 CondVector <- rep(timepoints,each=3)
 CondVector

# conditions must be specificed as factor, sorted along time-course 
Conditions <- factor(CondVector,levels=c("t3", "t7", "t10", "t13", "t15", "t17", "t19", "t21", "t23","t26"))
str(Conditions)
levels(Conditions)

######################

# move into directory where analysis will take place
setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/4_EBSeq_HMM/")


#####FLB analysis#####

# first need to get library size factors to adjust for seq depth differences between samples

Sizes <- MedianNorm(FLB_count_matrix)

# get the normalised matrix (so can plot graphs)
GeneNormData <- GetNormalizedMat(FLB_count_matrix,Sizes)

# now use EBSeqHMM to estimate gene expression and the probability that each gene is in a particular expression path (e.g. up-down-down-up)
# probably want to increase the UpdateRd parameter higher than 5 because will need more iterations

EBSeqHMMGeneOut_FLB <- EBSeqHMMTest(Data=FLB_count_matrix, sizeFactors=Sizes, Conditions=Conditions,UpdateRd=10)
# makes such a big file I don't think I would be able to use it!
saveRDS(EBSeqHMMGeneOut_FLB,"EBSeqHMMGeneOut_FLB.rds")

# detect DE genes under a target FDR (here FDR =0.001)
GeneDECalls_FLB <- GetDECalls(EBSeqHMMGeneOut_FLB, FDR=.001)
head(GeneDECalls_FLB)
# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_FLB,file="GeneDECalls_FLB_FDR0.001.csv")

# detect DE genes under a target FDR (here FDR =0.01)
GeneDECalls_FLB <- GetDECalls(EBSeqHMMGeneOut_FLB, FDR=.01)
head(GeneDECalls_FLB)
# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_FLB,file="GeneDECalls_FLB_FDR0.01.csv")

# detect DE genes under a target FDR (here FDR =0.05)
GeneDECalls_FLB <- GetDECalls(EBSeqHMMGeneOut_FLB, FDR=.05)
head(GeneDECalls_FLB)
# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_FLB,file="GeneDECalls_FLB_FDR0.05.csv")



# group gene paths together - keeping FDR =0.05 and cutoff for PP value 0.5, 
# might want to change OnlyDynamic to FALSE if want to also get genes which stay constant between some conditions
GeneConfCalls_FLB <- GetConfidentCalls(EBSeqHMMGeneOut_FLB, FDR=.05,cutoff=.5, OnlyDynamic=FALSE)

# write a CSV with only the confident called DE gene paths (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_FLB$Overall, file="GeneConfCalls_FLB_FDR0.05_cutoff0.5_only_dynamic_false.csv")

# write a CSV with only the number of confident called DE genes per path (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_FLB$NumEach, file="GeneConfCalls_FLB_num_each_FDR0.05_cutoff0.5_only_dynamic_false.csv")

## want to do some data visualisation

#plot QQ plot

pdf(file="QQP_FLB.pdf", height=10, width = 10)
par(mfrow=c(5,2))
QQP(EBSeqHMMGeneOut_FLB,GeneLevel=TRUE)
dev.off()

#plot DenNHist plot
pdf(file="DenNHist_FLB.pdf", height=10, width = 10)
par(mfrow=c(5,2))
DenNHist(EBSeqHMMGeneOut_FLB,GeneLevel=TRUE)
dev.off()


#####G analysis########
# first need to get library size factors to adjust for seq depth differences between samples

Sizes <- MedianNorm(G_count_matrix)

# get the normalised matrix (so can plot graphs)
GeneNormData <- GetNormalizedMat(G_count_matrix,Sizes)

# now use EBSeqHMM to estimate gene expression and the probability that each gene is in a particular expression path (e.g. up-down-down-up)
# probably want to increase the UpdateRd parameter higher than 5 because will need more iterations

EBSeqHMMGeneOut_G <- EBSeqHMMTest(Data=G_count_matrix, sizeFactors=Sizes, Conditions=Conditions,UpdateRd=10)
# makes such a big file I don't think I would be able to use it!
saveRDS(EBSeqHMMGeneOut_G,"EBSeqHMMGeneOut_G.rds")

# detect DE genes under a target FDR (here FDR =0.001)
GeneDECalls_G <- GetDECalls(EBSeqHMMGeneOut_G, FDR=.001)
head(GeneDECalls_G)
# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_G,file="GeneDECalls_G_FDR0.001.csv")

# detect DE genes under a target FDR (here FDR =0.01)
GeneDECalls_G <- GetDECalls(EBSeqHMMGeneOut_G, FDR=.01)
head(GeneDECalls_G)
# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_G,file="GeneDECalls_G_FDR0.01.csv")

# detect DE genes under a target FDR (here FDR =0.05)
GeneDECalls_G <- GetDECalls(EBSeqHMMGeneOut_G, FDR=.05)
head(GeneDECalls_G)
# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_G,file="GeneDECalls_G_FDR0.05.csv")

# group gene paths together - keeping FDR =0.05 and cutoff for PP value 0.5, 
# might want to change OnlyDynamic to FALSE if want to also get genes which stay constant between some conditions
GeneConfCalls_G <- GetConfidentCalls(EBSeqHMMGeneOut_G, FDR=.05,cutoff=.5, OnlyDynamic=FALSE)

# write a CSV with only the confident called DE gene paths (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_G$Overall, file="GeneConfCalls_G_FDR0.05_cutoff0.5_only_dynamic_false.csv")

# write a CSV with only the number of confident called DE genes per path (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_G$NumEach, file="GeneConfCalls_G_num_each_FDR0.05_cutoff0.5_only_dynamic_false.csv")

## want to do some data visualisation

#plot QQ plot

pdf(file="QQP_G.pdf", height=10, width = 10)
par(mfrow=c(5,2))
QQP(EBSeqHMMGeneOut_G,GeneLevel=TRUE)
dev.off()

#plot DenNHist plot
pdf(file="DenNHist_G.pdf", height=10, width = 10)
par(mfrow=c(5,2))
DenNHist(EBSeqHMMGeneOut_G,GeneLevel=TRUE)
dev.off()




