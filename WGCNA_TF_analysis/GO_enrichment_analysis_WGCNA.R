# Want to GO term enrichment for modules in WGCNA
# 26.9.2016
# Philippa Borrill

# read in module info

setwd("Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\maxP0.05\\")

module_info <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\maxP0.05\\genes_with_modules_mergeCutHeight0.15_no_expr_values.csv", 
                        header=T)
head(module_info)
dim(module_info)

# rename columns
colnames(module_info) <- c("gene", "bwnetModulelabels", "bwnetModuleColors")
head(module_info)
dim(module_info)

#rename rownames
rownames(module_info) <- module_info[,1]
module_info <- module_info[,-1]
head(module_info)

#### read in information about lengths and GO terms for each cluster #########

# read in GO terms
all_go <- read.csv("Y:\\PB_AFLF\\control_timecourse\\GO_terms_per_gene_long_form.csv")
head(all_go)
all_go <- all_go[,c(2,4)]
head(all_go)
dim(all_go)

# select only genes which were used for WGCNA
all_go <- subset(all_go, Gene %in% rownames(module_info))
dim(all_go)

#create vector for gene_lengths

# need to get lengths of genes not of transcripts
lengths <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\length_per_gene.csv", header=T)
head(lengths)
colnames(lengths) <- c("gene", "length")
head(lengths)

t1 <- subset(lengths, gene %in% rownames(module_info))
head(t1)
dim(t1)
# turn into a vector called gene.lens to use with GOSeq
gene.lens <- as.numeric(t1$length)
names(gene.lens) = t1$gene
head(gene.lens)


####### Do GO term enrichment and plot graph for each group ####
# start inspecting each cluster (this will save each hierachical group of genes and then perform GO term enrichment analysis using GOseq on each group
out_dir <- "Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\maxP0.05\\GO_enrichment_mergecutheight0.15\\"

assayed.genes <- as.vector(rownames(module_info))
length(assayed.genes)

library(goseq)
#i=1

for (i in seq(1, length(unique(module_info$bwnetModulelabels)))) {
  groupi_genes <- module_info[module_info$bwnetModulelabels == i, 1:2]
  groupi_genes$target_id <- rownames(groupi_genes)
  write.table(groupi_genes, file = paste0(out_dir, "module_", i, "_genes.tsv", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE)
  #get the GO terms
  ## Philippa thinks this step is unnecessary groupi_GO <- merge(groupi_genes, all_go, all.x = TRUE, all.y = FALSE)
  #now do GO stats analysis on the cluster using the cluster genes as the defined set and all genes used for clustering as the gene universe??
  #create a named binary vector for genes where one means expressed and 0 means not expressed
  de.genes <- groupi_genes$target_id
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  head(gene.vector)
  #now carry out the GOseq analysis
  pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
  GO.wall = goseq(pwf, gene2cat = all_go)
  #this gave table with p-values...now correct for multiple testing using FDR
  enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]
  head(enriched.GO)
  # add new column with over represented GO terms padj
  GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")
  # add new column with under represented GO terms padj
  GO.wall$under_rep_padj=p.adjust(GO.wall$under_represented_pvalue, method="BH")
  dep.GO=GO.wall$category[p.adjust(GO.wall$under_represented_pvalue, method="BH")<.05]
  write.table(GO.wall, file = paste0(out_dir, "module_", i, "_GOseq.tsv", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE)
  
}


# now extract GO over-representation info for revigo
setwd("Y:/PB_AFLF/control_timecourse/TF_analysis/WGCNA/maxP0.05/GO_enrichment_mergecutheight0.15")
getwd()

for (i in seq(1, length(unique(module_info$bwnetModulelabels)))) {
GO_data <- read.table(file = paste0(out_dir, "module_", i, "_GOseq.tsv", sep = ""), header=T, sep="\t")

head(GO_data)
dim(GO_data)
GO_data <- (GO_data[GO_data$over_rep_padj <0.05,c(1,8)])
dim(GO_data)
head(GO_data)

write.table(GO_data, file=paste0(out_dir, "module_", i, "_over_rep_GO_for_revigo.tsv", sep = ""), sep = "\t", quote=FALSE, col.names=F, row.names = F)
}

# now plot some graphs - heatmaps annotated with tissue 

install.packages("NMF")
library("NMF")


setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\maxP0.1\\heatmaps")

module_info2 <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\maxP0.1\\genes_with_modules_mergeCutHeight0.25.csv", header=T)

rownames(module_info2) <- module_info2$X
module_info2 <- module_info2[,-1]

# want to have some meta-data
metadata <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\metadata_for_WGCNA_expVIP_leaf.csv", header=T)
dim(metadata)
head(metadata)
dim(module_info2)
colnames(module_info2)
module_info2[1:4,1:10]
module_info2[1:4,326:340]
colnames(module_info2)
head(metadata)

i = 3

for (i in seq(3, length(unique(module_info2$bwnetModuleLabels)))) {


# select only orange module
module_sel <- module_info2[module_info2$bwnetModuleLabels == i,]
head(module_sel)

# get data as matrix
matrix_to_plot <- as.matrix(module_sel[,1:338])

# choose annotation for columns
annotation <- data.frame(metadata[,c(7:10)])
head(annotation)
#aheatmap(matrix_to_plot, annCol = annotation)

# choose annotation for rows
TF_family <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\triticum_aestivum_TFs_final_high_conf.csv", header=TRUE)
head(TF_family)
matrix_to_plot[1:4,1:4]

TF_annot <- merge(matrix_to_plot, TF_family, by.x=0, by.y="gene", all.x=T)
head(TF_annot)
dim(TF_annot)
TF_annot <- data.frame(TF_annot[,342])
head(TF_annot)
colnames(TF_annot) <- "family"
tail(TF_annot)


TF_annot$edited <- ifelse(TF_annot$family == "NAC", "NAC", 
                          ifelse(TF_annot$family == "WRKY", "WRKY", 
                                 ifelse( TF_annot$family == "MYB", "MYB",
                                         ifelse(TF_annot$family == "<NA>", "<NA>", "other_TF"))))
(TF_annot)

TF_annot <- data.frame(TF_annot$edited)
head(TF_annot)

nmf.options(grid.patch=TRUE)


pdf(file=paste0("module_",i,"_heatmap_ordered_tissue_age.pdf"), width = 12)
aheatmap(matrix_to_plot, annCol = annotation, annRow = TF_annot, Colv=order(metadata[,8],metadata[,9]) )
dev.off()

png(file=paste0("module_",i,"_heatmap_ordered_tissue_age.png"), width=4000, height =2000, res =250)
aheatmap(matrix_to_plot, annCol = annotation, annRow = TF_annot, Colv=order(metadata[,8],metadata[,9]))
dev.off()

# now plot with detailed tissue  info
# choose annotation for columns
annotation <- data.frame(metadata[,c(4)])
head(annotation)
colnames(annotation) <- "tissue"

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


pdf(file=paste0("module_",i,"_heatmap_detailed_ordered_tissue_age.pdf"), width = 12)
aheatmap(matrix_to_plot, annCol = annotation, annRow = TF_annot, annColors = list(col_vector), Colv=order(metadata[,8],metadata[,9]))
dev.off()

png(file=paste0("module_",i,"_heatmap_detailed_ordered_tissue_age.png"), width=4000, height =2000, res =250)
aheatmap(matrix_to_plot, annCol = annotation, annRow = TF_annot, annColors = list(col_vector), Colv=order(metadata[,8],metadata[,9]))
dev.off()


}

# what module are NAM and Y2H genes in?
getwd()
NAM_Y2H <- read.csv( "Y:/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/6_WGCNA_inc_expVIP_data/NAM_Y2H_genes.csv", header=F)
head(NAM_Y2H)

dim(NAM_Y2H)[1]

for (j in 1:dim(NAM_Y2H)[1]) {
print(  module_info2[as.character(NAM_Y2H[j,]),339:340] )
}
