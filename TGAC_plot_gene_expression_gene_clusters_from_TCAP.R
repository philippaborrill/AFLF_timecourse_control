# Aim to plot gene expression for particular genes
# 25.08.2016

# load meme results from TCAP which has group and gene info
#setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\7_TCAP_analysis\\low_and_high_conf\\functional_analysis_inputs")
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\7_TCAP_analysis\\Qian_-50\\functional_analysis_inputs")

meme <- read.table(file="meme.txt")
head(meme)
dim(meme)
colnames(meme) <- c("cluster", "gene")

# list of clusters to look at
cluster_list <- c("6", "15", "16", "19", "25", "86")

cluster_list <- seq(1,102)
cluster_list
cluster_list <- c("30","56")
cluster_list <- c("15")
cluster_list <- seq(1:43)
cluster_list <- seq(1:12)

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


#### can restart here to load vector of tpms per gene #####
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\1_SOM_analysis\\")

tpmData <- read.csv("tpm_summarised_per_gene.csv")
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

# now just select TF
TF_info <- read.csv(file="Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\gene_TF_family.csv")
head(TF_info)
colnames(TF_info) <- c("gene", "TF_family")
head(TF_info)

tpmF_TF <- merge(tpmF, TF_info, by.x ="target_id", by.y= "gene")
head(tpmF_TF)
dim(tpmF_TF)

# melt data to long form to use with ggplot
library(reshape2)

melted_tpmF <- melt(tpmF_TF, id.vars = c("target_id","TF_family"), variable.name = "sample")
head(melted_tpmF)
tail(melted_tpmF)
dim(melted_tpmF)

# now want to remove rep info from "sample" column so they can be easily averaged
melted_tpmF$timepoint <- gsub("_WT_.", "", (melted_tpmF$sample)) # . in the reg expression stands for any character
melted_tpmF$timepoint <- gsub("FLB", "", (melted_tpmF$timepoint))
head(melted_tpmF)
tail(melted_tpmF)



##### want to loop through each cluster and plot a graph for it #######

for (i in 1:length(cluster_list)) {
  setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\7_TCAP_analysis\\")
  print(i)
  j <- cluster_list[i]
  print(j)
  gene_list2 <- meme[which(meme$cluster==j),]
  head(gene_list2)
  # meme out put had TGACv1 with a capital V
  gene_list2$gene <- gsub("TGACV1", "TGACv1", (gene_list2$gene))
  head(gene_list2)


tpm_to_plot <- subset(melted_tpmF, melted_tpmF$target_id %in% gene_list2$gene)

head(tpm_to_plot)
dim(tpm_to_plot)
#write.csv(tpm_to_plot, paste0("cluster_",j,"_expression_TF_family.csv"))


# make summary of data


### summarise my data #######

summarised_tpm <- summarySE(tpm_to_plot, measurevar="value", groupvars=c("target_id","timepoint","TF_family"))

summarised_tpm

summarised_tpm$timepoint <- as.numeric(summarised_tpm$timepoint)

# setwd to save in

#setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\7_TCAP_analysis\\low_and_high_conf")
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\7_TCAP_analysis\\Qian_-10\\functional_analysis_inputs")


write.csv(summarised_tpm, paste0("cluster_",j,"_expression_TF_family.csv"))

# make graph
library(ggplot2)

plot <- ggplot(summarised_tpm, aes(x=timepoint, y=value, shape=TF_family, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point(size=3) +
  ylab("Expression (tpm)") +
  scale_x_continuous(breaks=seq(0,26,2)) +
  expand_limits(x=c(3,26)) + 
  guides(col = guide_legend(nrow = (2*(dim(genes)[1]))))
ggsave(plot, file = paste0("cluster_",j,"_all_genes.jpg"), width=20, height=20)
#ggsave(plot, file = paste0("cluster_",j,"_all_genes.jpg"), width=((dim(genes)[1])), height=((dim(genes)[1])),limitsize=F)

plot2 <- ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ylab("Expression (tpm)") +
  scale_x_continuous(breaks=seq(0,26,2)) +
  expand_limits(x=c(3,26)) +
  facet_wrap(~ TF_family, ncol=3 ) +
  ggtitle(paste0("cluster_",j,"_TCAP")) + 
  guides(col = guide_legend(nrow = (2*(dim(genes)[1]))))

ggsave(plot2, file = paste0("cluster_",j,"_per_TF_family.jpg"), width=20, height=((dim(genes)[1])))
#ggsave(plot2, file = paste0("cluster_",j,"_per_TF_family.jpg"), width=((dim(genes)[1])), height=((dim(genes)[1])),limitsize=F)

}
