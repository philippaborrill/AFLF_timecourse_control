# plot graph of EBSeqHMM normalised data for specific genes


# load packages
library("EBSeq")
library("EBSeqHMM")


# using the Vignette as a guide: http://bioconductor.org/packages/release/bioc/vignettes/EBSeqHMM/inst/doc/EBSeqHMM_vignette.pdf
# requires matrix of counts, with samples as columns, individual genes as rows
# the object "conditions" should be a factor with the same length as the number of samples e.g. 1 level per sample
# e.g. if you have 5 timepoints with 3 reps you would have t1,t1,t1,t2,t2,t2 etc
# I think I will need to analys the FLB and G samples separately

## read in data and only select genes with > 1 tpm ##

setwd("Y:/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/1_SOM_analysis/")

count.data<-read.csv(file="counts_summarised_per_gene.csv",header=T)  

#head(count.data)
dim(count.data)

# make rownames correct
rownames(count.data) <- count.data[,1]
counts <- count.data[,-1]
#head(counts)
print("dimensions of count data before filtering")
dim(counts)

# now only keep counts which are from TFs (exclude all other genes to speed up runtime)
# load file which only contains genes which are TFs
TF_family <- read.csv("Y:/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/gene_TF_family.csv", header =T)
head(TF_family)
colnames(TF_family) <- c("gene", "TF_family")
#head(TF_family)

# select only genes from tpmData_av which are TFs
counts <- merge(counts, TF_family, by.x = 0, by.y= "gene")
head(counts)
dim(counts)
dim(TF_family)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- counts[,1:60]
print("dimensions of count data TFs only")
dim(counts)
#head(counts)

# now select only genes with > 1tpm in at least 1 timepoint FLB

tpmData_FLB <- read.csv(file="tpm_summarised_per_gene.csv", header=T)
#head(tpmData_FLB)
colnames(tpmData_FLB)
# adjust rownames
rownames(tpmData_FLB) <- tpmData_FLB$X
tpmData_FLB <- tpmData_FLB[,-1]
#head(tpmData_FLB)

# just use FLB data
tpmData_FLB <- tpmData_FLB[,1:30]
#head(tpmData_FLB)


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
#head(tpmData_FLB)
colnames(tpmData_FLB)[31:40]
tpmData_FLB_av <- tpmData_FLB[,31:40]
#head(tpmData_FLB_av)

tpmData_FLB_av$maxtpm <- apply(tpmData_FLB_av[,1:10],1,max)
#head(tpmData_FLB_av)

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
#head(counts_FLB)
#head(row.names(counts_FLB))
# remove tpm and maxtpm columns
counts_FLB <- counts_FLB[,1:30]
#head(counts_FLB)
dim(counts_FLB)

#convert to matrix
count_matrix_FLB <- as.matrix(counts_FLB)
print("dimensions of count_matrix TFs in FLB >1 tpm:")
dim(count_matrix_FLB)

# rename to fit rest of script
FLB_count_matrix <- count_matrix_FLB


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
setwd("Y:/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/4b_EBSeq_HMM_TF_only/")


#####FLB analysis#####
library("EBSeqHMM")
# first need to get library size factors to adjust for seq depth differences between samples

Sizes <- MedianNorm(FLB_count_matrix)

# get the normalised matrix (so can plot graphs)
GeneNormData <- GetNormalizedMat(FLB_count_matrix,Sizes)

# plot graph 
plot.new()
par(mfrow=c(3,2))
PlotExp(GeneNormData, Conditions, Name="TRIAE_CS42_2AS_TGACv1_113243_AA0353410")
PlotExp(GeneNormData, Conditions, Name="TRIAE_CS42_2BS_TGACv1_145996_AA0452240")
PlotExp(GeneNormData, Conditions, Name="TRIAE_CS42_2DS_TGACv1_179582_AA0608070")
PlotExp(GeneNormData, Conditions, Name="TRIAE_CS42_6AS_TGACv1_486738_AA1564640")
PlotExp(GeneNormData, Conditions, Name="TRIAE_CS42_6DS_TGACv1_542626_AA1725630")


