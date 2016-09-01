# Aim to plot gene expression for particular genes
# 25.08.2016

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
head(melted_tpmF)
tail(melted_tpmF)

# now want to just select certain genes
# read in list of genes of interest
genes <- read.csv("2tpm\\2nd_run\\NAC_TF_in_group_5.csv", header=F)
genes
dim(genes)

# select only genes of interest in melted_tpmF
tpm_to_plot <- subset(melted_tpmF, target_id %in% genes$V1)
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

# make graph
library(ggplot2)

ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("pink", "hotpink", "violet", "#666666", "#66CC00", "#33CC66", "purple","#9966FF", "lightblue", "blue", "darkblue"))

#add the underlay lines
plot <- plot + geom_line(data = groupi.melted, mapping = aes(x = Time_point, y = value, group = variable), colour = "LIGHTGRAY", alpha = 1/2, size = 1/2)
#add the mean line
plot <- plot + geom_line(data = means.melted, mapping = aes(x = Time_point, y = value, group = variable, color = "Gene mean"), size=1)
plot <- plot + geom_line(data = hmeans.melted, mapping = aes(x = Time_point, y = value, group = variable, color = "Cluster mean") , size=1)
plot <- plot + scale_color_manual(values = c("Cluster mean" = "RED", "Gene mean" = "BLACK"))
plot <- plot + ggtitle(plot_title) + labs( y = "Normalised tpm") + annotate("text", label = paste0("n = ", n_transcripts), x= Inf , hjust = 1, y=Inf, vjust = 1)
ggsave(plot, file = plot_name)






##### make graphs for different gene sets ######
# read in list of genes of interest
genes <- read.csv("2tpm\\2nd_run\\Y2H_genes.csv", header=F)
genes
dim(genes)
head(genes)

# select only genes of interest in melted_tpmF
tpm_to_plot <- subset(melted_tpmF, target_id %in% genes$V1)
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

# make graph
library(ggplot2)

ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("pink", "hotpink", "violet", "#666666", "#66CC00", "#33CC66", "purple","#9966FF", "lightblue", "blue", "darkblue", "grey"), 
                       breaks = genes$V1, labels= paste(genes$V2,genes$V1)) +

plot <- ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("pink", "hotpink", "violet", "#666666", "#66CC00", "#33CC66", "purple","#9966FF", "lightblue", "blue", "darkblue", "grey"), 
                     breaks = genes$V1, labels= paste(genes$V2,genes$V1))
plot <- plot + ggtitle("Y2H_genes_tpm_expression")
plot <- plot + theme(legend.text = element_text(size=6))

ggsave(plot, file = "Y2H_genes_tpm.tiff",width=8, height=4)   




### replot without gene 5 ####

# read in list of genes of interest
genes <- read.csv("2tpm\\2nd_run\\Y2H_genes.csv", header=F)
genes
dim(genes)
genes <- genes[c(1:4,6:12),]
genes

# select only genes of interest in melted_tpmF
tpm_to_plot <- subset(melted_tpmF, target_id %in% genes$V1)
head(tpm_to_plot)
dim(tpm_to_plot)


# make summary of data

### summarise my data #######

summarised_tpm <- summarySE(tpm_to_plot, measurevar="value", groupvars=c("target_id","timepoint"))

summarised_tpm

# make graph
library(ggplot2)

ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("pink", "hotpink", "violet", "#666666", "#33CC66", "purple","#9966FF", "lightblue", "blue", "darkblue", "grey"), 
                     breaks = genes$V1, labels= paste(genes$V2,genes$V1)) 

plot <- ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("pink", "hotpink", "violet", "#666666",  "#33CC66", "purple","#9966FF", "lightblue", "blue", "darkblue", "grey"), 
                     breaks = genes$V1, labels= paste(genes$V2,genes$V1))
plot <- plot + ggtitle("Y2H_genes_tpm_expression")
plot <- plot + theme(legend.text = element_text(size=6))

ggsave(plot, file = "Y2H_genes_tpm_without_gene_5.tiff",width=8, height=4) 


### replot just gene 1+2 and 9+10 ####

# read in list of genes of interest
genes <- read.csv("2tpm\\2nd_run\\Y2H_genes.csv", header=F)
genes
dim(genes)
genes <- genes[c(1:2,8:9),]
genes

# select only genes of interest in melted_tpmF
tpm_to_plot <- subset(melted_tpmF, target_id %in% genes$V1)
head(tpm_to_plot)
dim(tpm_to_plot)


# make summary of data

### summarise my data #######

summarised_tpm <- summarySE(tpm_to_plot, measurevar="value", groupvars=c("target_id","timepoint"))

summarised_tpm

# make graph
library(ggplot2)

ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("pink", "hotpink", "#9966FF",  "blue"), 
                     breaks = genes$V1, labels= paste(genes$V2,genes$V1)) 

plot <- ggplot(summarised_tpm, aes(x=timepoint, y=value, colour=target_id, group =target_id)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("pink", "hotpink", "#9966FF",  "blue"), 
                     breaks = genes$V1, labels= paste(genes$V2,genes$V1))
plot <- plot + ggtitle("Y2H_genes 1, 2, 9 and 10")
plot <- plot + theme(legend.text = element_text(size=6))

ggsave(plot, file = "Y2H_genes_tpm_gene1_2_9_and_10.tiff",width=8, height=4) 

