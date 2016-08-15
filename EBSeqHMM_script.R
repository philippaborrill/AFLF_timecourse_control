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


setwd("/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/")

count.data<-read.table(file="final_output_counts.txt",header=T)  

head(count.data)
dim(count.data)

# make rownames correct
rownames(count.data) <- count.data[,1]
counts <- count.data[,-1]
head(counts)


#convert to matrix
count_matrix <- as.matrix(counts)
print("dimensions of count_matrix all genes:")
dim(count_matrix)

# get rid of genes with no expression in any sample
count_matrix <- count_matrix[(rowSums(count_matrix>0)>0),]
print("dimensions of count_matrix genes with zero expression removed:")
dim(count_matrix)

# separate FLB and G counts to two separate matrices
FLB_count_matrix <- count_matrix[,1:30]
head(FLB_count_matrix)
colnames(FLB_count_matrix)
dim(FLB_count_matrix)

G_count_matrix <- count_matrix[,31:60]
head(G_count_matrix)
colnames(G_count_matrix)
dim(G_count_matrix)

# make vector of time-points
timepoints <- c("t3", "t7", "t10", "t13", "t15", "t17", "t19", "t21", "t23","t26")
 CondVector <- rep(timepoints,each=3)
 CondVector

# conditions must be specificed as factor, sorted along time-course 
Conditions <- factor(CondVector,levels=c("t3", "t7", "t10", "t13", "t15", "t17", "t19", "t21", "t23","t26"))
str(Conditions)
levels(Conditions)

######################

#####FLB analysis#####

# first need to get library size factors to adjust for seq depth differences between samples

Sizes <- MedianNorm(FLB_count_matrix)

# get the normalised matrix (so can plot graphs)
GeneNormData <- GetNormalizedMat(FLB_count_matrix,Sizes)

# now use EBSeqHMM to estimate gene expression and the probability that each gene is in a particular expression path (e.g. up-down-down-up)
# probably want to increase the UpdateRd parameter higher than 5 because will need more iterations

EBSeqHMMGeneOut_FLB <- EBSeqHMMTest(Data=FLB_count_matrix, sizeFactors=Sizes, Conditions=Conditions,UpdateRd=10)
saveRDS(EBSeqHMMGeneOut_FLB,"EBSeqHMMGeneOut_FLB.rds")

# detect DE genes under a target FDR (here FDR =0.05)
GeneDECalls_FLB <- GetDECalls(EBSeqHMMGeneOut_FLB, FDR=.05)
head(GeneDECalls_FLB)

# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_FLB,file="GeneDECalls_FLB_FDR0.05.csv")

# group gene paths together - keeping FDR =0.05 and cutoff for PP value 0.5, 
# might want to change OnlyDynamic to FALSE if want to also get genes which stay constant between some conditions
GeneConfCalls_FLB <- GetConfidentCalls(EBSeqHMMGeneOut_FLB, FDR=.05,cutoff=.5, OnlyDynamic=TRUE)

# write a CSV with only the confident called DE gene paths (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_FLB$Overall, file="GeneConfCalls_FLB_FDR0.05_cutoff0.5_only_dynamic_true.csv")

# write a CSV with only the number of confident called DE genes per path (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_FLB$NumEach, file="GeneConfCalls_FLB_num_each_FDR0.05_cutoff0.5_only_dynamic_true.csv")

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
saveRDS(EBSeqHMMGeneOut_G,"EBSeqHMMGeneOut_G.rds")

# detect DE genes under a target FDR (here FDR =0.05)
GeneDECalls_G <- GetDECalls(EBSeqHMMGeneOut_G, FDR=.05)
head(GeneDECalls_G)

# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_G,file="GeneDECalls_G_FDR0.05.csv")

# group gene paths together - keeping FDR =0.05 and cutoff for PP value 0.5, 
# might want to change OnlyDynamic to FALSE if want to also get genes which stay constant between some conditions
GeneConfCalls_G <- GetConfidentCalls(EBSeqHMMGeneOut_G, FDR=.05,cutoff=.5, OnlyDynamic=TRUE)

# write a CSV with only the confident called DE gene paths (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_G$Overall, file="GeneConfCalls_G_FDR0.05_cutoff0.5_only_dynamic_true.csv")

# write a CSV with only the number of confident called DE genes per path (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_G$NumEach, file="GeneConfCalls_G_num_each_FDR0.05_cutoff0.5_only_dynamic_true.csv")

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




