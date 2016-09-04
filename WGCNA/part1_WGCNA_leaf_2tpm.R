# aim is to run WGCNA on leaf samples from RNA-seq
# 31-08-2016
# Philippa Borrill

#Using the tutorials at https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/ as  guide

# first install the packages which I need
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute","AnnotationDbia")) 
install.packages("WGCNA") 


##### Tutorial 1 ########
# Using tutorial 1 as a guide for this section
# In the tutorial they are using microarray data so I need to do a bit more pre-processing to my data

# Peter Langfelder wrote that you should: 
# 1) "filter out genes whose count is less than say 5 in more than say 80 % of samples" this gets rid of a lot of noise 
# I have now decide to filter by tpm because even filtering by counts left many low expressed genes
# 2) "use a variance-stabilizing transformation, such as the one implemented in varianceStabilizingTransformation or rlogTransformation in DESeq2

# Therefore I Want to use count data which is required by DESeq and will normalise using DESeq 

#### Loading data and pre-processing ####
#set working directory
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\2_Leaf_WGCNA\\2_tpm_1_timepoint\\")

# load in count data which has already been summarised per gene by tximport from the script "self_organising_map_analysis.R"
count.data <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\counts_summarised_per_gene.csv")
head(count.data)
dim(count.data)

# this includes all the data including the grain but I just want to look at leaf samples for now
# therefore just keep leaf sample columns

leaf.count.data <- count.data[,1:31]
dim(leaf.count.data)

# make rownames correct
rownames(leaf.count.data) <- leaf.count.data[,1]
counts <- leaf.count.data[,-1]
head(counts)
head(row.names(counts))

# round counts to integers (required for DESeq)
counts <- round(counts)
dim(counts)

# Filter data to only keep high confidence genes with over certain tpm expression in at least 1 sample
# tpm filter level
tpm <- 3

# First select high confidence genes therefore need to get list of high confidence genes
gene_conf <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_confidence.csv", header=TRUE)
head(gene_conf)
colnames(gene_conf) <- c("gene","confidence")
head(gene_conf)
# add confidence level to counts
counts_conf <- merge(counts, gene_conf, by.x = 0, by.y= "gene")
head(counts_conf)
dim(counts_conf)
# select only high confidence genes
counts_conf <- counts_conf[which(counts_conf$confidence == "High"),]
# remove column saying confidence level
counts_conf <- counts_conf[,1:31]
head(counts_conf)
dim(counts_conf)


# clean up workspace to remove unnecessary dataframes
rm(count.data)
rm(counts)
rm(gene_conf)
rm(leaf.count.data)

# now select only genes with > 2tpm in at least 1 timepoint
#### can restart here to load vector of tpms per gene #####
tpmData <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\tpm_summarised_per_gene.csv")
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

# clean up workspace to remove unnecessary dataframes
rm(tpmData)

# merge together counts_conf  and tpmData_av
counts_conf_max <- merge(counts_conf, tpmData_av, by.x = "Row.names", by.y = 0)
head(counts_conf_max)
dim(counts_conf_max)
# select only rows with a maxtpm >2
counts <- counts_conf_max[which(counts_conf_max$maxtpm>tpm),]

# make rownames correct
rownames(counts) <- counts[,1]
counts <- counts[,-1]
head(counts)
head(row.names(counts))
# remove tpm and maxtpm columns
counts <- counts[,1:30]
head(counts)
dim(counts)

# now prepare for DESeq matrix loading
# assign levels to columns to use as factors
timepoints <- c("FLB3", "FLB7", "FLB10", "FLB13", "FLB15", "FLB17", "FLB19", "FLB21", "FLB23","FLB26")
CondVector <- rep(timepoints,each=3)
CondVector

samples <- data.frame(row.names=colnames(counts), condition=as.factor(CondVector))
samples

library(DESeq2)
# get data into DESeq form
dds <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design=~condition)
head(dds)

#make sure the FLB3 is used as the reference condition :
colData(dds)$condition <- relevel(colData(dds)$condition, "FLB3")

# run Deseq
dds2 <- DESeq(dds)
head(dds2)

#do variance stabilising transformation (see https://www.biostars.org/p/95788/)
vsd_blind <- varianceStabilizingTransformation(dds2,blind=TRUE)
vc <- assay(vsd_blind)
head(vc)

# transform the data to be samples as rows and genes as columns (what WGCNA expects)
vsd2 <- t(vc)
vsd2[1:4,1:4]
dim(vsd2)


######## starting WGCNA analysis ########
# source the WGCNA package
library("WGCNA")
options(stringsAsFactors = FALSE)
# if running in R-studio need to disable multi-threading
disableWGCNAThreads()

# if running in R-gui can use multi-threading
allowWGCNAThreads()

# copy vsd2 to a new dataframe called dataExpr0 to match the manual
datExpr0 <- vsd2

# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK 

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
#if (!gsg$allOK)
#   {if (sum(!gsg$goodGenes)>0)
#       printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
#       if (sum(!gsg$goodSamples)>0)
#           printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
#       datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
#       }

#All my genes were OK

# Now check all samples look ok 
sampleTree = hclust(dist(datExpr0), method = "average");

pdf(file = paste0("Sample_Clustering_",tpm,"tpm.pdf"), width = 6, height = 4)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# If want to add trait data (e.g. chlorophyll level or age or protein content or moisture content need to do this now)
# see 1c in https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf 

# save the data to use in the next step
save(datExpr0, file=paste0("filtered_leaf_data_ready_for_WGCNA_",tpm,"tpm.RData"))
