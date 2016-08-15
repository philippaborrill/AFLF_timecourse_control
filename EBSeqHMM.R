# Script to run EBSeqHMM to look at genes differentially expressed across time in RNA-seq timecourse
# Philippa Borrill 08.02.2016



# install packages
source("https://bioconductor.org/biocLite.R")
biocLite("EBSeq")
biocLite("EBSeqHMM")
biocLite("DESeq")

# load packages
library("EBSeq")
library("EBSeqHMM")


# using the Vignette as a guide: http://bioconductor.org/packages/release/bioc/vignettes/EBSeqHMM/inst/doc/EBSeqHMM_vignette.pdf
# requires matrix of counts, with samples as columns, individual genes as rows
# the object "conditions" should be a factor with the same length as the number of samples e.g. 1 level per sample
# e.g. if you have 5 timepoints with 3 reps you would have t1,t1,t1,t2,t2,t2 etc
# I think I will need to analys the FLB and G samples separately

## if using cluster
source R-3.1.0
R
setwd("/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/")

# if using desktop
setwd("Y:\\PB_AFLF\\control_timecourse\\kallisto_analysis\\add_samples")


count.data<-read.table(file="final_output_counts.txt",header=T)  

head(count.data)
dim(count.data)

# make rownames correct
rownames(count.data) <- count.data[,1]
counts <- count.data[,-1]
head(counts)

#convert to matrix
count_matrix <- as.matrix(counts)
head(count_matrix)
colnames(count_matrix)

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

# plot expression of 2 genes:
PlotExp(GeneNormData, Conditions, Name="Traes_7DS_FFA36F6DA.1")
PlotExp(GeneNormData, Conditions, Name="Traes_7DS_FFE9ACDAB.2")

# now use EBSeqHMM to estimate gene expression and the probability that each gene is in a particular expression path (e.g. up-down-down-up)
# probably want to increase the UpdateRd parameter higher than 5 because will need more iterations

EBSeqHMMGeneOut_FLB <- EBSeqHMMTest(Data=FLB_count_matrix, sizeFactors=Sizes, Conditions=Conditions,UpdateRd=5)

# ran out of memory on desktop for EBSeqHMMGeneOut function therefore will need to use the cluster - might try interactive node


#####G analysis########
# first need to get library size factors to adjust for seq depth differences between samples

Sizes_G <- MedianNorm(G_count_matrix)

# get the normalised matrix (so can plot graphs)
GeneNormData_G <- GetNormalizedMat(G_count_matrix,Sizes_G)

PlotExp(GeneNormData_G, Conditions, Name="Traes_7DS_FFA36F6DA.1")
PlotExp(GeneNormData_G, Conditions, Name="Traes_7DS_FFE9ACDAB.2")




##########
## make practise set with fewer genes
FLB_count_matrix_100 <- FLB_count_matrix[2000:2100,]
head(FLB_count_matrix_100)

# first need to get library size factors to adjust for seq depth differences between samples

Sizes <- MedianNorm(FLB_count_matrix_100)

# get the normalised matrix (so can plot graphs)
GeneNormData <- GetNormalizedMat(FLB_count_matrix_100,Sizes)

# plot expression of 2 genes:
PlotExp(GeneNormData, Conditions, Name="Traes_1AL_4338FCC28.2")

# now use EBSeqHMM to estimate gene expression and the probability that each gene is in a particular expression path (e.g. up-down-down-up)
# probably want to increase the UpdateRd parameter higher than 5 because will need more iterations

EBSeqHMMGeneOut_FLB <- EBSeqHMMTest(Data=FLB_count_matrix_100, sizeFactors=Sizes, Conditions=Conditions,UpdateRd=5)

# detect DE genes under a target FDR (here FDR =0.05)
GeneDECalls_FLB <- GetDECalls(EBSeqHMMGeneOut_FLB, FDR=.05)
head(GeneDECalls_FLB)

# write csv of differentially expressed gene paths (one row per gene)
write.csv(GeneDECalls_FLB,file="GeneDECalls_FLB.csv")

# group gene paths together - keeping FDR =0.05 and cutoff for PP value 0.5, 
# might want to change OnlyDynamic to FALSE if want to also get genes which stay constant between some conditions
GeneConfCalls_FLB <- GetConfidentCalls(EBSeqHMMGeneOut_FLB, FDR=.05,cutoff=.5, OnlyDynamic=TRUE)

# write a CSV with only the confident called DE gene paths (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(GeneConfCalls_FLB$Overall, file="GeneConfCalls_FLB_FDR0.5_cutoff0.5_only_dynamic_true.csv")

# write a CSV with only the number of confident called DE genes per path (e.g. filtered by FDR, cutoff PP and dynamic path)
write.csv(Up-Down-Up-Down-Up-Up-Up-Down-Down$NumEach, file="GeneConfCalls_FLB_num_each_FDR0.5_cutoff0.5_only_dynamic_true.csv")

## want to do some data visualisation
GeneNormData <- GetNormalizedMat(FLB_count_matrix_100, Sizes)
print(GeneConfCalls_FLB$EachPath[["Up-Down-Up-Down-Up-Up-Up-Down-Down"]])

GeneOfInterest <- GeneConfCalls_FLB$EachPathNames[["Up-Down-Up-Down-Up-Up-Up-Down-Down"]]
print(GeneOfInterest)

par(mfrow=c(3,1))
for(i in 1:3)PlotExp(FLB_count_matrix_100, Conditions, Name=GeneOfInterest[i])

pdf(file="QQP_FLB.pdf", height=10, width = 10)
par(mfrow=c(5,2))
QQP(EBSeqHMMGeneOut_FLB,GeneLevel=TRUE)
dev.off()


pdf(file="DenNHist_FLB.pdf", height=10, width = 10)
par(mfrow=c(5,2))
DenNHist(EBSeqHMMGeneOut_FLB,GeneLevel=TRUE)
dev.off()

