# aim is to run WGCNA on leaf samples from RNA-seq
# 31-08-2016
# Philippa Borrill

#Using the tutorials at https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/ as  guide

#### TUTORIAL 5 ########
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.pdf 
# network visualisation

# setup R correctly
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\2_Leaf_WGCNA\\")
library(WGCNA)
options(stringsAsFactors = FALSE)

#load data from the 1st part of the analysis
lnames=load(file="filtered_leaf_data_ready_for_WGCNA.RData")
# lnames contains the names of loaded variables
lnames

datExpr <- datExpr0
rm(datExpr0)

# Load network data saved in the second part.
# let's use the network from:deepSplit =4, detectCutHeight = 0.99995, mergeCutHeight =0.25 
# This is the one where the NAM genes aren't all in the same cluster

lnames=load(file="detectCutHeight_0.99995\\bwnet2_network_components_mergeCutHeight0.25_detectCutHeight_0.99995_deepSplit_4.RData")
lnames

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


## first part is to plot the TOM per block ##
# takes a long time so skip for now

MEs = moduleEigengenes(datExpr, bwnetModuleColors)$eigengenes
head(MEs)

# order MEs (at this point could add in trait data)
MET = orderMEs(MEs)

# plot the dendrogram

jpeg(file="Eigengene_dendrogram_detectCutHeight_0.99995_mergeCutHeight_0.25_deepSplit_4.jpg", height=600, width=800)
par(cex=1)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro=c(0,4,2,0),
                      plotHeatmaps=FALSE)
dev.off()

jpeg(file="Eigengene_adjacency_heatmap_detectCutHeight_0.99995_mergeCutHeight_0.25_deepSplit_4.jpg", height=600, width =600)
par(cex=1)
plotEigengeneNetworks(MET,"Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms=FALSE, xLabelsAngle=90)
dev.off()


######## Tutorial 6########## Simulated data (III)
# assessing quality of modules #

datME = moduleEigengenes(datExpr, bwnetModuleColors)$eigengenes
signif(cor(datME, use="p"), 2)

dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram

jpeg(file="Eigengene_dendrogram_detectCutHeight_0.99995_mergeCutHeight_0.25_deepSplit_4_tutorial-6.jpg", height=600, width=800)

par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")
dev.off()

# pairswise scatter plots of the samples (arrays) along the module eigengenes
plot.new()
jpeg(file="Scatter_plot_eigengenesdetectCutHeight_0.99995_mergeCutHeight_0.25_deepSplit_4_tutorial-6.jpg", height=6000, width=8000)
par(oma=c(1,1,1,1))
par(mar=c(1,1,1,1))
plotMEpairs(datME)
dev.off()
