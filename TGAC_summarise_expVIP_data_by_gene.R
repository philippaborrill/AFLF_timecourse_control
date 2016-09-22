# Aim: to summarise transcript expression per gene for expVIP samples selected for WGCNA

wd <- "Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data\\"

read_dir <- "Y:\\expression_browser\\TGAC_assembly\\analysis\\kallisto_results_bootstrap\\results"

source("https://bioconductor.org/biocLite.R")
biocLite("tximportData")
install.packages("readr")
library(tximportData)
library(readr)

# make vector pointing to the kallisto results files   ########
samples <- read.table(paste0(wd,"samples.txt"), header=T)
samples

files <- file.path(read_dir, samples$sample, "abundance.tsv", fsep ="\\")
files
names(files) <- paste0(samples$sample)
files
all(file.exists(files))

# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.table("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\transcripts_to_genes.txt", header=T)
head(tx2gene)

library(tximport)
# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
names(txi)

# move into directory where I will save this analysis
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\6_WGCNA_inc_expVIP_data")

# to see counts summarised per gene
head(txi$counts)
colnames(txi$counts)


# save counts summarised per gene
write.csv(txi$counts, file="counts_summarised_per_gene_expVIP.csv")

# to see tpm summarised per gene
head(txi$abundance)
colnames(txi$abundance)

# save tpm summarised per gene
write.csv(txi$abundance, file="tpm_summarised_per_gene_expVIP.csv")

# see lengths summarised per gene
head(txi$length)

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="length_per_gene_expVIP.csv")
