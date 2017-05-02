# aim is to run WGCNA on expVIP samples with TGAC genesq
# 31-08-2016 # edited 20-04-2017
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
setwd("Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\")

# read expVIP data
counts <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\counts_summarised_per_gene_expVIP.csv")
head(counts)
dim(counts)


# make rownames correct
rownames(counts) <- counts[,1]
counts <- counts[,-1]
head(counts)
head(row.names(counts))

# round counts to integers (required for DESeq)
counts <- round(counts)
dim(counts)


# Filter data to only keep high confidence genes with over certain tpm expression in at least 1 sample
# tpm filter level
tpm_threshold <- 0.5

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
counts_conf <- counts_conf[,1:309]
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


# now select only genes with > 0.5tpm in at least 1 sample
#### can restart here to load vector of tpms per gene #####

# now add expVIP data
tpm.data.expVIP <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\tpm_summarised_per_gene_expVIP.csv")
head(tpm.data.expVIP)
dim(tpm.data.expVIP)

# make rownames correct
rownames(tpm.data.expVIP) <- tpm.data.expVIP[,1]
tpm <- tpm.data.expVIP[,-1]
head(tpm)
head(row.names(tpm))
dim(tpm)

library(matrixStats)

tpm_filt <- tpm[rowCounts(as.matrix(tpm>tpm_threshold))>=3,] # select only rows which have expr >0.5 tpm in >=3 samples
head(tpm)
head(tpm_filt)
dim(tpm_filt)


# clean up workspace to remove unnecessary dataframes
rm(tpm.data.expVIP)


# merge together counts_conf  and tpm_filt # when merging only rows with >0.5 tpm in 3 samples will be kept
counts_conf_max <- merge(counts_conf, tpm_filt, by.x = 0, by.y = 0)
head(counts_conf_max)
dim(counts_conf_max)
counts <- counts_conf_max

# make rownames correct
rownames(counts) <- counts[,1]
counts <- counts[,-1]
head(counts)
head(row.names(counts))
# remove tpm  columns
dim(counts)
colnames(counts[1:308])
counts <- counts[,1:308]
colnames(counts) <-gsub(".x","",colnames(counts))

head(counts)
dim(counts)

class(counts)

#clean up workspace
rm(counts_conf_max)
rm(counts_conf)
rm(tpm)
rm(tpm_filt)


source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("lattice")
biocLite("colorspace")

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

pdf(file = paste0("Sample_Clustering_",tpm_threshold,"tpm.pdf"), width = 20, height = 4)
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
metadata_selected <- metadata_selected[,-3]
head(metadata_selected)
dim(metadata_selected)


# save the data to use in the next step
save(datExpr0, metadata_selected, file=paste0("filtered_expVIP_data_ready_for_WGCNA_high_conf_",tpm_threshold,"tpm.RData"))



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

pdf(file = paste0("Sample_Clustering_coloured_by_tissue_",tpm_threshold,"tpm.pdf"), width = 40, height = 4)
par(cex = 0.6);
par(mar = c(6,4,2,0))
plot(dend, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


# now plot with high level age as colour
metadata_selected$High.level.age
colorCodes <- c(seedling = "pale green", vegetative = "green3", reproductive= "olivedrab")
labels_colors(dend) <- colorCodes[(metadata_selected$High.level.age[(order.dendrogram(dend))])]
labels_colors(dend)
head(metadata_selected$High.level.age)
head(metadata_selected$High.level.age[(order.dendrogram(dend))])

labels_colors(dend)

pdf(file = paste0("Sample_Clustering_coloured_by_age_",tpm_threshold,"tpm.pdf"), width = 40, height = 4)
par(cex = 0.6);
par(mar = c(6,4,2,0))
plot(dend, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


