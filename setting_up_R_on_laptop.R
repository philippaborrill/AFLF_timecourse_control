# Aim is to set up and install packages required for SOM and WGCNA analysis
# 04.09.2016
# Philippa Borrill

setwd("C:/Users/borrillp/Documents/AFLF_control_timecourse/")
getwd()

# find out where libraries will be saved
.libPaths()
#[1] "C:/Users/borrillp/Documents/R/win-library/3.3"
#[2] "C:/Program Files/R/R-3.3.1/library"

#install tximport for importing from kallisto
source("https://bioconductor.org/biocLite.R")
biocLite("tximportData")
biocLite("tximport")
install.packages("readr")

# install EBSeqHMM
source("https://bioconductor.org/biocLite.R")
biocLite("EBSeqHMM")

# install timecourse
source("https://bioconductor.org/biocLite.R")
biocLite("timecourse")
browseVignettes("timecourse")

# install DESeq2
biocLite("DESeq2")

# install SOM
install.packages("kohonen")

# install WGCNA
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute","AnnotationDbia")) 
install.packages("WGCNA") 

install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("ggplot2")
biocLite("goseq")



