# Want to plot eigengenes for modules in WGCNA
# 05.05.2017
# Philippa Borrill

# check we can load the data
setwd("Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\maxP0.05\\")

library(WGCNA)
options(stringsAsFactors = FALSE)

#load data from the 1st part of the analysis
lnames=load(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\filtered_expVIP_data_ready_for_WGCNA_high_conf_0.5tpm.Rdata")
# lnames contains the names of loaded variables
lnames

datExpr <- datExpr0
rm(datExpr0)

# Load network data saved in the second part.
# mergecutheight0.15

lnames=load(file="bwnet_network_components_mergeCutHeight0.15.Rdata")
lnames

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

MEs = moduleEigengenes(datExpr, bwnetModuleLabels)$eigengenes

# or if want colours rather than numbers:
#MEs = moduleEigengenes(datExpr, bwnetModuleColors)$eigengenes

head(MEs)
write.csv(MEs, file="eigengenes.csv")

#plot deodrogram of Eigengene relatedness
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# default order of MEs is:
rownames(datExpr)

# need to re-order into something sensible to be able to spot patterns  and perhaps rename samples (e)

names(MEs)

#reorder according to a sensible order I worked out in excel arrange first by tissue and then at the end it has the stress samples
# saved as "eigengenes_with_metadata.xsxl" tab "reordered_metadata"
MEs_reordered <- MEs[c(135,136,81,137,138,165,169,170,175,176,181,188,190,191,192,193,194,199,200,
                       205,206,211,218,219,222,223,228,229,234,235,240,241,248,166,167,173,174,177,
                       178,182,183,186,187,195,196,201,202,207,208,212,213,220,221,224,225,230,231,
                       236,237,242,243,139,140,168,171,172,179,180,184,185,189,197,198,203,204,209,
                       210,214,215,216,217,226,227,232,233,238,239,244,245,246,247,85,86,87,113,116,
                       120,50,52,112,114,117,54,17,36,38,39,40,42,14,16,28,35,18,25,27,29,51,53,13,15,
                       19,30,115,118,119,55,21,22,26,31,41,32,33,34,20,23,24,37,82,83,84,153,154,155,
                       156,157,158,88,89,90,100,101,102,121,122,141,142,57,58,59,60,7,8,9,143,144,161,
                       162,65,66,67,68,73,74,75,76,159,160,163,164,145,146,47,48,49,43,44,45,46,147,
                       148,1,2,3,149,150,69,70,71,72,77,78,79,80,151,152,103,104,105,106,107,108,109,
                       110,111,91,92,93,94,95,96,97,98,99,56,63,64,61,62,10,11,12,4,5,6,123,124,131,
                       132,127,128,125,126,133,134,129,130,250,252,259,263,267,269,271,273,275,283,
                       293,295,297,299,305,251,256,257,258,261,262,276,277,280,282,284,288,298,300,
                       307,253,254,260,264,265,278,279,281,286,289,292,296,302,306,308,249,255,266,
                       268,270,272,274,285,287,290,291,294,301,303,304),]

MEs[304,]

MEs_reordered[308,]
# reordering worked as expected
names_order <- read.csv(file="names_and_order_for_MEs.csv", header=T)
head(names_order)

rownames(MEs_reordered)

MEs_reordered$names <- names_order$Name.for.R

barplot(MEs_reordered$MEred)
barplot(MEs_reordered$MEtan, names.arg = MEs_reordered$names)

rownames(MEs)

plot(MEs_reordered$MEred)
lines(MEs_reordered$MEred)
axis(1, labels=T)

# try melting data to plot with ggplot
library(reshape2)
melted_MEs <- melt(MEs_reordered)
head(melted_MEs)

colnames(melted_MEs) <- c("tissue", "module", "value")
head(melted_MEs)

# need to make meltedMEs$tissue an ordered factor
melted_MEs$tissue <- factor(melted_MEs$tissue, levels= melted_MEs$tissue)

library(ggplot2)

pdf(file="module_eigengene_expr.pdf", width=10, height=80)

ggplot(data=melted_MEs, aes(x=tissue, y=value)) +
  geom_bar(stat="identity") + 
    facet_wrap(~ module, ncol=1) + theme_minimal() + theme(axis.line= element_line(color="black"), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#plot heatmap of MEs
head(MEs_reordered)
rownames(MEs_reordered) <- MEs_reordered$names
dim(MEs_reordered)

pdf(file="module_eigengene_expr_heatmap.pdf", width =20, height=10)
heatmap(t(as.matrix(MEs_reordered[,1:38])),Colv =NA, cexCol = 0.1)
dev.off()
