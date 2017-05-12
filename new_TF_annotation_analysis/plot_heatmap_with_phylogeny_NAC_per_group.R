# Aim is to plot dendrogram of NAC phylogeny with expression data next to it as a heatmap
# 09-05-2017
# Philippa Borrill


#### Loading data and pre-processing ####
#set working directory

# decide to use data which was already filtered for WGCNA: i.e. only includes genes with >0.5 tpm in at least 3 samples
# the count data was used and this was normalised using varianceStabilizingTransformation from DESeq2
setwd("Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\maxP0.05\\")

install.packages("Matrix")
install.packages("foreign")
install.packages("rpart")
install.packages("htmltools")
library(WGCNA)
options(stringsAsFactors = FALSE)

#load data from the 1st part of the analysis
lnames=load(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\filtered_expVIP_data_ready_for_WGCNA_high_conf_0.5tpm.Rdata")
# lnames contains the names of loaded variables
lnames

head(metadata_selected)

datExpr <- datExpr0
rm(datExpr0)

datExpr[1:4,1:4]

# transform datExpr to have samples as columns and genes as columns
t_datExpr <- t(datExpr)
t_datExpr[1:4,1:4]
rownames(t_datExpr)
# make into a dataframe
t_datExpr <- as.data.frame(t_datExpr)
t_datExpr[1:4,1:4]

t_datExpr_reordered <- t_datExpr[,c(135,136,81,137,138,165,169,170,175,176,181,188,190,191,192,193,194,199,200,
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
                                        268,270,272,274,285,287,290,291,294,301,303,304)]





### if want to reload count data and re-normalise
# load in count data which has already been summarised per gene by tximport from the script "self_organising_map_analysis.R"
count.data <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\counts_summarised_per_gene_expVIP.csv")
head(count.data)
dim(count.data)

# make rownames correct
rownames(count.data) <- count.data[,1]
counts <- count.data[,-1]
head(counts)
head(row.names(counts))
is.data.frame(counts)

# want to normalise so expression is on a scale of 0 to 1 for each gene
counts_norm <- t(apply(counts, 1, function(x)(x-min(x))/(max(x)-min(x))))
head(counts_norm)
counts_norm <- as.data.frame(counts_norm)
dim(counts_norm)

# want to re-order the columns of counts_norm to be something sensible:
counts_norm_reordered <- counts_norm[,c(135,136,81,137,138,165,169,170,175,176,181,188,190,191,192,193,194,199,200,
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
                       268,270,272,274,285,287,290,291,294,301,303,304)]


#try log(counts)

log_counts <- log(counts) 

log_counts_reordered <- log_counts[,c(135,136,81,137,138,165,169,170,175,176,181,188,190,191,192,193,194,199,200,
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
                                        268,270,272,274,285,287,290,291,294,301,303,304)]


# try tpm
# load in count data which has already been summarised per gene by tximport from the script "self_organising_map_analysis.R"
tpm.data <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\tpm_summarised_per_gene_expVIP.csv")
head(tpm.data)
dim(tpm.data)

# make rownames correct
rownames(tpm.data) <- tpm.data[,1]
tpms <- tpm.data[,-1]
head(tpms)
head(row.names(tpms))
is.data.frame(tpms)

library(matrixStats)
tpm_filt <- tpms[rowCounts(as.matrix(tpms>0.5))>=3,] # select only rows which have expr >0.5 tpm in >=3 samples

# want to normalise so expression is on a scale of 0 to 1 for each gene
tpms_filt_norm <- t(apply(tpm_filt, 1, function(x)(x-min(x))/(max(x)-min(x))))
head(tpms_filt_norm)
tpms_filt_norm <- as.data.frame(tpms_filt_norm)
dim(tpms_filt_norm)

# want to re-order the columns of tpms_norm to be something sensible:
tpms_norm_reordered <- tpms_filt_norm[,c(135,136,81,137,138,165,169,170,175,176,181,188,190,191,192,193,194,199,200,
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
                                        268,270,272,274,285,287,290,291,294,301,303,304)]




# load in tree
install.packages("nlme")
library("ape")
library("Biostrings")
library("ggplot2")
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("ggtree")
library("ggtree")
biocLite("ggrepel")
library("ggrepel")

# move into working directory
setwd("Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\phylogenetics\\trees\\groups_a-h")
getwd()

# load in tree file
tree <- read.tree(file="NAC-a_wheat_longID.phy.tree")

tree
plot(tree)

# need to edit tree to have gene ID not transcript ID
tree$tip.label <- gsub("\\.\\d","",tree$tip.label)

list_genes <- as.data.frame(tree$tip.label)
list_genes
colnames(list_genes) <- c("gene")
head(list_genes)


ggtree(tree) + geom_tiplab()
# add a scale
ggtree(tree) + geom_treescale()

# plot as cladogram
ggtree(tree, branch.length="none")

#### try adding in heatmap info ###
p <- ggtree(tree, branch.length="none") + geom_tiplab()
p

p <- ggtree(tree, branch.length="none") + geom_tiplab(cex = 1) + xlim(0,150)
print(p)

# if need to check column names:
#p2 <- gheatmap(p, counts_norm_reordered,  low="yellow", high="red",offset=15, color=NULL, width=5, colnames_position="bottom",  font.size=2 ) +xlim(0,100)

#column names removed because hard to read
# using varianceStabilizingTransformation from DESeq2:
p2 <- gheatmap(p, t_datExpr_reordered,  low="yellow", high="red",offset=15, color=NULL, width=5, colnames=F ) +xlim(0,100)
print(p2)

pdf(file="test_datExpr_NAC_phy_full_geneIDs_heatmap_norm_counts.pdf", width = 10, height =5)
print(p2)
dev.off()

# using normalisation across a gene to 1:
p3 <- gheatmap(p, counts_norm_reordered,  low="yellow", high="red",offset=15, color=NULL, width=5, colnames=F ) +xlim(0,100)
print(p3)

pdf(file="test_counts_norm_NAC_phy_full_geneIDs_heatmap_norm_counts.pdf", width = 10, height =5)
print(p3)
dev.off()

# using log counts:
p4 <- gheatmap(p, log_counts_reordered,  low="yellow", high="red",offset=15, color=NULL, width=5, colnames=F ) +xlim(0,100)
print(p4)

pdf(file="test_log_counts_NAC_phy_full_geneIDs_heatmap_norm_counts.pdf", width = 10, height =5)
print(p4)
dev.off()


####### now want to plot graphs for each group a-h ###########
# load in tree file
groups <- c("a","b","c","d","e","f","g","h")

for (i in groups){
  print(i)
  tree <- read.tree(file=paste0("NAC-",i,"_wheat_longID.phy.tree"))
  
  # need to edit tree to have gene ID not transcript ID
  tree$tip.label <- gsub("\\.\\d","",tree$tip.label)
  
  p <- ggtree(tree, branch.length="none") + geom_tiplab(cex = 1) + xlim(0,150)
  
  # using varianceStabilizingTransformation from DESeq2:
  p2 <- gheatmap(p, t_datExpr_reordered,  low="yellow", high="red",offset=15, color=NULL, width=5, colnames=F ) +xlim(0,100)
  print(p2)
  
  pdf(file=paste0("VSD_NAC-",i,"_heatmap.pdf"), width = 10, height =5)
  print(p2)
  dev.off()
  
  # using normalisation across a gene to 1:
  p3 <- gheatmap(p, counts_norm_reordered,  low="yellow", high="red",offset=15, color=NULL, width=5, colnames=F ) +xlim(0,100)
  print(p3)
  
  pdf(file=paste0("counts_norm_NAC-",i,"_heatmap.pdf"), width = 10, height =5)
  print(p3)
  dev.off()
  
  
  }


for (i in groups){
  print(i)
  tree <- read.tree(file=paste0("NAC-",i,"_wheat_longID.phy.tree"))
  
  # need to edit tree to have gene ID not transcript ID
  tree$tip.label <- gsub("\\.\\d","",tree$tip.label)
  
  p <- ggtree(tree, branch.length="none") + geom_tiplab(cex = 1) + xlim(0,150)
  
  # using log2(tpm):
  p2 <- gheatmap(p, tpms_norm_reordered,  low="yellow", high="red",offset=15, color=NULL, width=5, colnames=F ) +xlim(0,100)
  print(p2)
  
  pdf(file=paste0("tpm_filt_low_NAC-",i,"_heatmap.pdf"), width = 10, height =5)
  print(p2)
  dev.off()
  
  
  
  
}


