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
# 2) "use a variance-stabilizing transformation, such as the one implemented in varianceStabilizingTransformation or rlogTransformation in DESeq2

# Therefore I Want to use count data which is required by DESeq and will normalise using DESeq 


#### Loading data and pre-processing ####
#set working directory
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\2_Leaf_WGCNA")

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

#Peter Langfelder suggested to filter out genes whose count is less than say 5 in more than say 80 % of samples

# make new data frame which has the same format as counts but is a logical of 1 count is over 5 or 0 if count is under 5
logical.counts <- counts>5
head(counts)
head(logical.counts)
# get rowsums of logical.counts (i.e. how many samples have expr over 5 counts)
num_samples_over_5 <- apply(logical.counts,1,sum)
head(num_samples_over_5)

# now just want to keep genes which have over 5 expression in at least 20 % of samples (6 samples)
counts.filt20 <- counts[num_samples_over_5>5,]
head(counts.filt20)
dim(counts.filt20)

# now just want to keep genes which have over 5 expression in at least 10 % of samples (3 samples)
counts.filt10 <- counts[num_samples_over_5>2,]
head(counts.filt10)
dim(counts.filt10)
rm(counts.filt10)

# now prepare for DESeq matrix loading
# assign levels to columns to use as factors
timepoints <- c("FLB3", "FLB7", "FLB10", "FLB13", "FLB15", "FLB17", "FLB19", "FLB21", "FLB23","FLB26")
CondVector <- rep(timepoints,each=3)
CondVector

samples <- data.frame(row.names=colnames(counts.filt20), condition=as.factor(CondVector))
samples

library(DESeq2)
# get data into DESeq form
dds <- DESeqDataSetFromMatrix(countData = counts.filt20, colData=samples, design=~condition)
head(dds)

#make sure the FLB3 is used as the reference condition :
colData(dds)$condition <- relevel(colData(dds)$condition, "FLB3")

# run Deseq
dds2 <- DESeq(dds)
head()

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

pdf(file = "Sample_Clustering.pdf", width = 6, height = 4)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# If want to add trait data (e.g. chlorophyll level or age or protein content or moisture content need to do this now)
# see 1c in https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf 

# save the data to use in the next step
save(datExpr0, file="filtered_leaf_data_ready_for_WGCNA.RData")
