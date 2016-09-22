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
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\")

# load in count data which has already been summarised per gene by tximport from the script "self_organising_map_analysis.R"
count.data <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\counts_summarised_per_gene.csv")
head(count.data)
dim(count.data)

# this includes all the data including the grain but I just want to look at leaf samples for now
# therefore just keep leaf sample columns
leaf.count.data <- count.data[,1:31]
dim(leaf.count.data)

# now add expVIP data
count.data.expVIP <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\counts_summarised_per_gene_expVIP.csv")
head(count.data.expVIP)
dim(count.data.expVIP)

# merge expVIP and leaf data
merged_counts <- merge(count.data.expVIP, leaf.count.data, by.x= "X", by.y="X" )
head(merged_counts)
dim(merged_counts)

# make rownames correct
rownames(merged_counts) <- merged_counts[,1]
counts <- merged_counts[,-1]
head(counts)
head(row.names(counts))

# round counts to integers (required for DESeq)
counts <- round(counts)
dim(counts)

# remove unnecessary datasets
rm(leaf.count.data)
rm(count.data.expVIP)
rm(count.data)

# Filter data to only keep high confidence genes with over certain tpm expression in at least 1 sample
# tpm filter level
tpm <- 2

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
counts_conf <- counts_conf[,1:339]
head(counts_conf)
dim(counts_conf)

# make rownames correct
rownames(counts_conf) <- counts_conf[,1]
counts_conf <- counts_conf[,-1]
head(counts_conf)
head(rownames(counts_conf))
dim(counts_conf)

# clean up workspace to remove unnecessary dataframes
rm(counts)
rm(gene_conf)


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

# now add expVIP data
tpm.data.expVIP <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\tpm_summarised_per_gene_expVIP.csv")
head(tpm.data.expVIP)
dim(tpm.data.expVIP)

# merge expVIP and leaf data
merged_tpm <- merge(tpm.data.expVIP, tpmData_av, by.x= "X", by.y=0 )
head(merged_tpm)
dim(merged_tpm)

# make rownames correct
rownames(merged_tpm) <- merged_tpm[,1]
tpm <- merged_tpm[,-1]
head(tpm)
head(row.names(tpm))
dim(tpm)

tpmData_av$maxtpm <- apply(tpm[,1:318],1,max)
head(tpmData_av)

# clean up workspace to remove unnecessary dataframes
rm(tpmData)
rm(tpm)

# merge together counts_conf  and tpmData_av
counts_conf_max <- merge(counts_conf, tpmData_av, by.x = 0, by.y = 0)
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
dim(counts)
colnames(counts[1:338])
counts <- counts[,1:338]
head(counts)
dim(counts)

class(counts)




library(DESeq2)
# get data into DESeq form (convert to matrix)
counts_matrix <- as.matrix(counts)
head(counts_matrix)

#do variance stabilising transformation (see https://www.biostars.org/p/95788/)
vsd_blind <- varianceStabilizingTransformation(counts_matrix,blind=TRUE)
vc <- (vsd_blind)
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

pdf(file = paste0("Sample_Clustering_",tpm,"tpm.pdf"), width = 20, height = 4)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# If want to add trait data (e.g. chlorophyll level or age or protein content or moisture content need to do this now)
# should add in TISSUE + AGE info
# see 1c in https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf 
metadata <- read.csv(file="metadata_for_wGCNA.csv")
head(metadata)
dim(metadata)

metadata_selected <- subset(metadata, run_accession %in% rownames(datExpr0))
head(metadata_selected)
dim(metadata_selected)
rownames(metadata_selected) <- metadata_selected$run_accession
metadata_selected <- metadata_selected[,-2]
head(metadata_selected)
dim(metadata_selected)


# save the data to use in the next step
save(datExpr0, metadata_selected, file=paste0("filtered_leaf_data_ready_for_WGCNA_",tpm,"tpm.RData"))


##### couldn't get this to work! #####
# now plot dendrogram with colours according to tissue
install.packages("dendextend")
library("dendextend")

metadata_selected$High.level.tissue <- gsub("leaves/shoots", "leavesshoots", metadata_selected$High.level.tissue)
(metadata_selected$High.level.tissue)

# rename 

# convert sampleTree to dendrogram



dend <- as.dendrogram(sampleTree)



colorCodes <- c(roots= "brown", leavesshoots = "green", spike = "orange", grain= "purple")
labels_colors(dend) <- colorCodes[(metadata_selected$High.level.tissue[(order.dendrogram(dend))])]
labels_colors(dend)
head(metadata_selected$High.level.tissue)
head(metadata_selected$High.level.tissue[(order.dendrogram(dend))])

labels_colors(dend)

pdf(file = paste0("Sample_Clustering_coloured_by_tissue_",tpm,"tpm.pdf"), width = 40, height = 4)
par(cex = 0.6);
par(mar = c(6,4,2,0))
plot(dend, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()




