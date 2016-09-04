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
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\2_Leaf_WGCNA\\3_tpm_1_timepoint")
library(WGCNA)
options(stringsAsFactors = FALSE)

#load data from the 1st part of the analysis
lnames=load(file="filtered_leaf_data_ready_for_WGCNA_3tpm.RData")
# lnames contains the names of loaded variables
lnames
#check datExpr0 looks ok:
datExpr0[1:4,1:4]
dim(datExpr0)
# rename datExpr0 to match the tutorial
datExpr <- datExpr0
rm(datExpr0)
datExpr[1:4,1:4]
dim(datExpr)

# realise I should have run the pickSoftThreshold function with a signed network with bicor correlation therefore re-run
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# run soft-thresholding
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType ="signed hybrid", 
                        corFnc ="bicor", corOptions=list(maxPOutliers=0.1))


# Plot the results:
sizeGrWindow(9, 5)
pdf(file="soft-threshold_power_signed_hybrid_3tpm.pdf")
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

# This shows that for a signed hybrid network I should use a power of 9. 
power <- 9

##Run the model

bwnet = blockwiseModules(datExpr, maxBlockSize = 10000,
                         power = power, networkType = "signed hybrid", TOMType = "unsigned", minModuleSize = 30,
                         corType="bicor", corOptions = "use = 'p', maxPOutliers = 0.1",
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Leaf_power9_signed_hybrid_TOM-blockwise",
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


# check we can load the data
bwnet_test <- load(file="bwnet_network_mergeCutHeight0.25.RData")
bwnet_test


# looks like my cutoff for module membership was way too low so many genes remain unassigned 
# to get more modules need to reduce the merge cut height e.g. to 0.1


bwnet2 <- recutBlockwiseTrees(datExpr, goodSamples=bwnet$goodSamples, goodGenes =bwnet$goodGenes,
                    blocks = bwnet$blocks, TOMFiles = bwnet$TOMFiles, dendrograms = bwnet$dendrograms,
                    corType = "bicor", corOptions = "use = 'p', maxPOutliers = 0.1", networkType = "signed hybrid",
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
write.csv(gene_expr,"genes_with_modules_mergeCutHeight0.25.csv")

#write.csv(gene_expr,"genes_with_modules_mergeCutHeight0.25_detectCutHeight_0.9995_deepSplit_4.csv")

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

