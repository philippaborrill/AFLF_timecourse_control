# aim is to plot bar charts to show TF distribution between co-expression clusters

# Philippa Borrill 10-05-2017

setwd("Y:\\PB_AFLF\\control_timecourse\\TF_analysis\\WGCNA\\maxP0.05\\GO_enrichment_mergecutheight0.15\\")

mod_data <- read.csv(file="gene_modules_TFs_for_R.csv", header =T)
head(mod_data)

gsub("-","NA",mod_data)

library(plyr)

s_mod_data <- ddply(mod_data, c("Module"), summarise,
               "Number"    = length(Gene),
)
head(s_mod_data)



library(ggplot2)
totals <- ggplot(s_mod_data, aes(x= Module, y= Number))  + geom_bar(stat="identity")
totals

TF_mod_data <- mod_data[mod_data$TF.family != "-",]
head(TF_mod_data)
dim(TF_mod_data)

library(reshape2)
melted_mod_data <- melt(mod_data, id.vars = "Module",)
head(melted_mod_data)





TFs <- ggplot(TF_mod_data, aes(x = factor(Module), y= Gene, color = TF.family)) + geom_bar(stat ="identity")
TFs

ggplot(TF_mod_data)
