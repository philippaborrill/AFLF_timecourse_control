# aim is to run WGCNA on leaf samples from RNA-seq
# 31-08-2016
# Philippa Borrill

#Using the tutorials at https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/ as  guide

#### TUTORIAL 2c ########
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.pdf 
# Automatic network construction and module detection with block-wise network construction
# My data set had >10,000 probes which is approx the max that an 8Gb memory computer can handle
# Therefore I need to split hte network into blocks to make it possible to compute
# I guess I could run it on the cluster if I wanted

# setup R correctly
# if in cluster
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\")

library(WGCNA)
options(stringsAsFactors = FALSE)

##########network analysis itself ##########
#load data from the 1st part of the analysis
lnames=load(file="filtered_leaf_data_ready_for_WGCNA_2tpm.RData")
# lnames contains the names of loaded variables
lnames
#check datExpr0 looks ok:
datExpr0[1:4,1:4]

# rename datExpr0 to match the tutorial
datExpr <- datExpr0
rm(datExpr0)
datExpr[1:4,1:4]
dim(datExpr)

# change into new folder to save results
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\maxP0.1")


# realise I should have run the pickSoftThreshold function with a signed network with bicor correlation therefore re-run
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# run soft-thresholding
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType ="signed hybrid", 
                        corFnc ="bicor", corOptions=list(maxPOutliers=0.1))


# Plot the results:
sizeGrWindow(9, 5)
pdf(file="soft-threshold_power_signed_hybrid_2tpm.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# This shows that for a signed hybrid network I should use a power of 10. 
power <- 10


##Run the model
# cannot allocate vector to memory!
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000, 
                         power = power, networkType = "signed hybrid", TOMType = "unsigned", minModuleSize = 30,
                         corType="bicor", corOptions = "use = 'p', maxPOutliers = 0.1", # change to 0.09 so can't exclude a whole timepoint
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Leaf_power9_signed_hybrid_TOM-blockwise_maxP0.1",
                         verbose = 3)

# Examine module membership

# First let's see what the bwnet dataset contains
names(bwnet)
length(table(bwnet$blocks))
table(bwnet$blocks)

# get the modules colours
bwnetModuleColors <- labels2colors(bwnet$colors)

#get the module labels
bwnetModuleLabels <- bwnet$colors

# look at how many genes per module
table(bwnetModuleLabels)
write.csv(table(bwnetModuleLabels), file="bwnet_modules_mergeCutHeight_0.25.csv")

# plot histogram of module sizes
jpeg(file="Number of genes in module mergeCutHeight_0.25.jpg")
par(mfrow=c(1,2))
hist(table(bwnetModuleLabels), xlim=c(1,10000), breaks=1000, xlab="Number of genes in Module", main="mergeCutHeight_0.25")
hist(table(bwnetModuleLabels), xlim=c(1,1000), breaks=1000, xlab="Number of genes in Module", main="Zoomed in mergeCutHeight_0.25")

dev.off()

# plot dendrograms of module membership

# as 1 pdf (1 page per block)
pdf(file="dendrogram_of_module_membership_0.25_mergeCutHeight.pdf", width=6, height=4)

for (i in 1:length(table(bwnet$blocks)))  {        # there are 5 blocks
plotDendroAndColors(bwnet$dendrograms[[i]], bwnetModuleColors[bwnet$blockGenes[[i]]],
                    "Module colors", main = paste0("Gene dendrogram and module colors in block ",i),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
}
dev.off()


# as 1 jpg per block
for (i in 1:length(table(bwnet$blocks)))  {        # there are 5 blocks
  jpeg(file=paste0("dendrogram_of_module_membership_0.25_mergeCutHeight_block_",i,".jpeg"))
  plotDendroAndColors(bwnet$dendrograms[[i]], bwnetModuleColors[bwnet$blockGenes[[i]]],
                      "Module colors", main = paste0("Gene dendrogram and module colors in block ",i),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

# let's save this data
save(bwnet, file="bwnet_network_mergeCutHeight0.25.RData")


# also save the components separately (as the manual suggests)
bwnetMEs <- bwnet$MEs
bwnetdendrograms <- bwnet$dendrograms

save(bwnetModuleColors, bwnetModuleLabels, bwnetMEs, bwnetdendrograms, file = "bwnet_network_components_mergeCutHeight0.25.RData")

######playing with clustering settings ########
# check we can load the data
setwd("C:\\Users\\borrillp\\Documents\\AFLF_control_timecourse\\2_tpm_1_timepoint\\maxP0.09\\")
bwnet_test <- load(file="bwnet_network_mergeCutHeight0.25.RData")
bwnet_test


# looks like my cutoff for module membership was way too low so many genes remain unassigned 
# to get more modules need to reduce the merge cut height e.g. to 0.1

# need to rename TOM files because original calculation was with a different name
bwnet$TOMFiles
bwnet$TOMFiles2 <- as.vector(c("Leaf_power9_signed_hybrid_TOM-blockwise-block.1.RData",
                               "Leaf_power9_signed_hybrid_TOM-blockwise-block.2.RData",
                               "Leaf_power9_signed_hybrid_TOM-blockwise-block.3.RData",
                               "Leaf_power9_signed_hybrid_TOM-blockwise-block.4.RData",
                               "Leaf_power9_signed_hybrid_TOM-blockwise-block.5.RData"))
bwnet$TOMFiles <- bwnet$TOMFiles2

bwnet2 <- recutBlockwiseTrees(datExpr, goodSamples=bwnet$goodSamples, goodGenes =bwnet$goodGenes,
                    blocks = bwnet$blocks, TOMFiles = bwnet$TOMFiles, dendrograms = bwnet$dendrograms,
                    corType = "bicor", corOptions = "use = 'p', maxPOutliers = 0.09", networkType = "signed hybrid",
                    minModuleSize = 30, reassignThreshold = 0, 
                    mergeCutHeight = 0.25, # have changed this from 0.25 to 0.1 
                    detectCutHeight = 0.9995, # changed from the default of 0.995
                    deepSplit = 4 , # changed from default of 2
                    numericLabels = TRUE,
                    verbose = 3)
                    
# First let's see what the bwnet dataset contains
names(bwnet2)

# get the modules colours
bwnet2ModuleColors <- labels2colors(bwnet2$colors)

#get the module labels
bwnet2ModuleLabels <- bwnet2$colors

# look at how many genes per module
table(bwnet2ModuleLabels)
write.csv(table(bwnet2ModuleLabels), file="bwnet2_modules_mergeCutHeight_0.2mergeCutHeight_0.25_detectCutHeight_0.9995_deepSplit_4.csv")

# plot histogram of module sizes
jpeg(file="Number of genes in module mergeCutHeight_0.25_detectCutHeight_0.9995_deepSplit_4.jpg")
par(mfrow=c(1,2))
hist(table(bwnet2ModuleLabels), xlim=c(1,10000), breaks=1000, xlab="Number of genes in Module", main="mergeCutHeight_0.25
     _detectCutHeight_0.9995
     _deepSplit_4")
hist(table(bwnet2ModuleLabels), xlim=c(1,1000), breaks=1000, xlab="Number of genes in Module", main="Zoomed in mergeCutHeight_0.25_
     detectCutHeight_0.9995
     _deepSplit_4")

dev.off()

# plot dendrograms of module membership

# as 1 pdf (10 pages)
pdf(file="dendrogram_of_module_membership_0.25_mergeCutHeight_detectCutHeight_0.9995_deepSplit_4.pdf", width=6, height=4)
par(mfrow = c(5,2))
for (i in 1:10)  {        # there are 10 blocks
  plotDendroAndColors(bwnet$dendrograms[[i]], bwnet2ModuleColors[bwnet$blockGenes[[i]]],
                      "Module colors", main = paste0("Gene dendrogram and module colors in block ",i),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}
dev.off()


# as 10 jpgs
for (i in 1:10)  {        # there are 10 blocks
  jpeg(file=paste0("dendrogram_of_module_membership_0.25_mergeCutHeight_detectCutHeight_0.9995_deepSplit_4_block_",i,".jpeg"))
  plotDendroAndColors(bwnet$dendrograms[[i]], bwnet2ModuleColors[bwnet$blockGenes[[i]]],
                      "Module colors", main = paste0("Gene dendrogram and module colors in block ",i),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}                

# now let's save bwnet2 #NB that doesn't include the dendrogram- need that from bwnet original file
save(bwnet2, file="bwnet2_network_mergeCutHeight0.25_detectCutHeight_0.9995_deepSplit_4.RData")
names(bwnet2)

# also save the components separately (as the manual suggests)
bwnetMEs <- bwnet2$MEs
bwnetdendrograms <- bwnet$dendrograms
# get the modules colours
bwnetModuleColors <- labels2colors(bwnet2$colors)

#get the module labels
bwnetModuleLabels <- bwnet2$colors

# 
save(bwnetModuleColors, bwnetModuleLabels, bwnetMEs, bwnetdendrograms, file = 
       "bwnet2_network_components_mergeCutHeight0.25_detectCutHeight_0.9995_deepSplit_4.RData")


### now want to output genes with module membership
gene_expr <- cbind(t(datExpr),bwnetModuleLabels, bwnetModuleColors)
head(gene_expr)

# write genes with modules to file
write.csv(gene_expr,"genes_with_modules_mergeCutHeight0.25_detectCutHeight_0.9995_deepSplit_4.csv")
write.csv(gene_expr,"genes_with_modules_mergeCutHeight0.25.csv")


# what modules are the NAM genes in?

gene_expr["TRIAE_CS42_6AS_TGACv1_486738_AA1564640",31:32]
gene_expr["TRIAE_CS42_6BS_TGACv1_513229_AA1635270",31:32]
gene_expr["TRIAE_CS42_6DS_TGACv1_542626_AA1725630",31:32]
gene_expr["TRIAE_CS42_2AS_TGACv1_113243_AA0353410",31:32]
gene_expr["TRIAE_CS42_2BS_TGACv1_145996_AA0452240",31:32]
gene_expr["TRIAE_CS42_2DS_TGACv1_179582_AA0608070",31:32]

# all except 2A are unassigned suggesting they have low variance - plot to take a look

pdf(file="NAM_gene_expression_vsd_transformed.pdf",width=12, height=6)
par(mfrow=c(2,3))
par(oma=c(2,1,1,1))
names_plot <- colnames(gene_expr)[1:30]
names_plot
barplot(as.numeric(gene_expr["TRIAE_CS42_6AS_TGACv1_486738_AA1564640",1:30]), main="6A")
barplot(as.numeric(gene_expr["TRIAE_CS42_6BS_TGACv1_513229_AA1635270",1:30]), main="6B")
barplot(as.numeric(gene_expr["TRIAE_CS42_6DS_TGACv1_542626_AA1725630",1:30]), main="6D")
barplot(as.numeric(gene_expr["TRIAE_CS42_2AS_TGACv1_113243_AA0353410",1:30]), main="2A", names.arg = names_plot, las=2)
barplot(as.numeric(gene_expr["TRIAE_CS42_2BS_TGACv1_145996_AA0452240",1:30]), main="2B", names.arg = names_plot, las=2)
barplot(as.numeric(gene_expr["TRIAE_CS42_2DS_TGACv1_179582_AA0608070",1:30]), main="2D", names.arg = names_plot, las=2)

dev.off()

# now want to see what transcription factors are in the same modules as the NAM genes

setwd("C:\\Users\\borrillp\\Documents\\AFLF_control_timecourse\\2_tpm_1_timepoint\\maxP0.09\\")
loaded_modules <- read.csv(file="genes_with_modules_mergeCutHeight0.25.csv")
head(loaded_modules)
rownames(loaded_modules) <- loaded_modules$X
loaded_modules <- loaded_modules[,-1]
head(loaded_modules)
dim(loaded_modules)

# add in TF information
TF_info <- read.csv(file="C:\\Users\\borrillp\\Documents\\AFLF_control_timecourse\\gene_TF_family.csv")
head(TF_info)
colnames(TF_info) <- c("gene", "TF_family")
head(TF_info)

# merge TF info with the module info
merge_genes <- merge(loaded_modules, TF_info, by.x = 0, by.y = "gene",all.x=TRUE) # keep all genes in loaded_modules
head(merge_genes)
dim(merge_genes)

# save modules with TF info (not info about expression)
merge_genes_sel_col <- merge_genes[,c(1,32:34)]
head(merge_genes_sel_col)
write.csv(merge_genes_sel_col,file="genes_with_module_and_TF_info.csv")

# count number of genes in sienna3 module
length(merge_genes_sel_col$bwnetModuleColors[which(merge_genes_sel_col$bwnetModuleColors=="sienna3")])

# select only genes in module with NAM 6A (and 6D and 2D) (sienna 3)
NAM_6A <- merge_genes_sel_col[which(merge_genes_sel_col$bwnetModuleColors=="sienna3"),]
head(NAM_6A)
dim(NAM_6A)

# select only TFs
NAM_6A_TF <- NAM_6A[which(NAM_6A$TF_family != "<NA>"),]
NAM_6A_TF
dim(NAM_6A_TF)
write.csv(NAM_6A_TF, file="genes_in_sienna3_module35_which_are_TFs.csv")



#### now plot expression of TF genes in sienna3 #####
#### can restart here to load vector of tpms per gene #####

tpmData <- read.csv("C:\\Users\\borrillp\\Documents\\AFLF_control_timecourse\\tpm_summarised_per_gene.csv")
head(tpmData)
colnames(tpmData)
# adjust rownames
rownames(tpmData) <- tpmData$X
tpmData <- tpmData[,-1]
head(tpmData)

# just use FLB data
tpmData <- tpmData[,1:30]
head(tpmData)

# add rownames as new column of gene names
tpmData$target_id <- rownames(tpmData)
head(tpmData)
colnames(tpmData)
tpmF <- tpmData[, c(31, 1:30)]
rownames(tpmF) <- NULL
dim(tpmF)
head(tpmF)

# melt data to long form to use with ggplot
library(reshape2)

melted_tpmF <- melt(tpmF, id.vars = "target_id", variable.name = "sample")
head(melted_tpmF)
tail(melted_tpmF)
dim(melted_tpmF)

# now want to remove rep info from "sample" column so they can be easily averaged
melted_tpmF$timepoint <- gsub("_WT_.", "", (melted_tpmF$sample)) # . in the reg expression stands for any character
head(melted_tpmF)
tail(melted_tpmF)

# now want to just select certain genes
# read in list of genes of interest
genes <- read.csv("genes_in_sienna3_module35_which_are_TFs.csv", header=T)
genes
dim(genes)

# select only genes of interest in melted_tpmF
tpm_to_plot <- subset(melted_tpmF, target_id %in% genes$Row.names)
head(tpm_to_plot)
dim(tpm_to_plot)


# make summary of data

#### function to summarise data ##########
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
### summarise my data #######

summarised_tpm <- summarySE(tpm_to_plot, measurevar="value", groupvars=c("target_id","timepoint"))

summarised_tpm

# set name to save plot 
plot_name <- "modules35_sienna3_TF_gene_expression_in_tpm.pdf"

# make graph
library(ggplot2)
pdf(plot_name, height = 6, width=12)

ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("pink", "hotpink", "violet", "#666666", "#66CC00", "#33CC66", "purple","#9966FF", "lightblue", "blue", "darkblue"))
dev.off()


