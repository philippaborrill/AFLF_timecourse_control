# Aim to plot gene expression for particular genes
# 25.08.2016

# gene list to plot is csv 
# don't have header
genes_to_plot <- "NAC_TF_in_group_5"
genes_to_plot <- "NAM_genes"
genes_to_plot <- "NAC_in_cluster_30_TCAP"
genes_to_plot <- "genes_for_crossing_sept_2016"
genes_to_plot <- "genes_to_cross_11"
genes_to_plot <- "genes_with_3_or_more_connections"

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

# melt data to long form to use with ggplot
library(reshape2)

melted_tpmF <- melt(tpmF, id.vars = "target_id", variable.name = "sample")
head(melted_tpmF)
tail(melted_tpmF)
dim(melted_tpmF)

# now want to remove rep info from "sample" column so they can be easily averaged
melted_tpmF$timepoint <- gsub("_WT_.", "", (melted_tpmF$sample)) # . in the reg expression stands for any character
melted_tpmF$timepoint <- gsub("FLB", "", (melted_tpmF$timepoint))
head(melted_tpmF)
tail(melted_tpmF)

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


  
  wd <- "Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\5_plot_expression\\"
  gene_list <- paste0(wd,genes_to_plot,".csv")
  
  

# now want to just select certain genes
# read in list of genes of interest
genes <- read.csv(gene_list, header=F)
genes
dim(genes)
dim(genes)[1]

# select only genes of interest in melted_tpmF
tpm_to_plot <- subset(melted_tpmF, target_id %in% genes$V1)
head(tpm_to_plot)
dim(tpm_to_plot)
### summarise my data #######

summarised_tpm <- summarySE(tpm_to_plot, measurevar="value", groupvars=c("target_id","timepoint"))

summarised_tpm

summarised_tpm$timepoint <- as.numeric(summarised_tpm$timepoint)

# setwd to save in


# make graph
library(ggplot2)

ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  ylab("Expression (tpm)") +
  scale_x_continuous(breaks=seq(0,26,2)) +
  expand_limits(x=c(3,26))
  
plot <- ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  ylab("Expression (tpm)") +
  scale_x_continuous(breaks=seq(0,26,2)) +
  expand_limits(x=c(3,26))
ggsave(plot, file = paste0(wd,genes_to_plot,"_all_genes.jpg"), width=10, height=6)

plot2 <- ggplot(summarised_tpm, aes(x=timepoint, y=value)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ylab("Expression (tpm)") +
  scale_x_continuous(breaks=seq(0,26,2)) +
  expand_limits(x=c(3,26)) +
  facet_wrap(~ target_id, ncol=3 ) +
  ggtitle(genes_to_plot)

ggsave(plot2, file = paste0(wd,genes_to_plot,"_per_genes.jpg"), width=10, height=dim(genes)[1])


