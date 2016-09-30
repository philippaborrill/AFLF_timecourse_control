# script to select genes for CSI for 10 tpm in specific families

# Philippa Borrill
# 30.09.2016


setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\9_EBSeqHMM_TF_only_new_annotation\\")

# read file which was produced by the scipt TGAC_EBSeqHMM_TF_only_try_different_tpm_and_prep_data_TCAP_CSI.R
tpm_data <- read.csv("countData_diff_expr_TF_10tpm_high_conf_with_TF.csv", header=T)
head(tpm_data)

# now read in file containing list of NAC, MYB and WRKY which were found in CSI network 10 tpm 0.014 threshold
TF_in_network <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\8_CSI_analysis\\new_annotation_NAC_WRKY_MYB\\10_tpm\\0.014\\0.014_node_table.csv", header=T)
head(TF_in_network)


###### Extract data for two batch approach: ######

setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\8_CSI_analysis\\new_annotation_NAC_WRKY_MYB\\all_families")

two_batch_run1 <- read.csv("two_batch_run1.csv", header=F)
head(two_batch_run1)
two_batch_run2 <- read.csv("two_batch_run2.csv", header=F)
head(two_batch_run2)

# select only genes which are in the TF_network from CSI
data_TF_in_network <- merge(tpm_data, TF_in_network, by.x = "condition", by.y = "name")
head(data_TF_in_network)
colnames(data_TF_in_network)[1:31]
data_TF_in_network <- data_TF_in_network[,1:31]
head(data_TF_in_network)
dim(data_TF_in_network)

# select only genes which are in the TF families selected for two_batch_run1
colnames(tpm_data)

run1 <- tpm_data[which(tpm_data$family %in% two_batch_run1$V1),1:31]
head(run1)
dim(run1)

# add together the NAC,MYB,WKRY data with the run1 data
run1_complete <- rbind.data.frame(data_TF_in_network, run1)
head(run1_complete)
dim(run1_complete)

# add in correct column names
colnames(run1_complete) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(run1_complete)
run1_complete <- rbind(colnames(run1_complete), run1_complete)
colnames(run1_complete) <- c("condition",rep("WT",30))
head(run1_complete)

# write csv
write.csv(run1_complete, file="two_batches_run1_for_CSI.csv",row.names = F)


# select only genes which are in the TF families selected for two_batch_run2
run2 <- tpm_data[which(tpm_data$family %in% two_batch_run2$V1),1:31]
head(run2)
dim(run2)

# add together the NAC,MYB,WKRY data with the run2 data
run2_complete <- rbind.data.frame(data_TF_in_network, run2)
head(run2_complete)
dim(run2_complete)

# add in correct column names
colnames(run2_complete) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(run2_complete)
run2_complete <- rbind(colnames(run2_complete), run2_complete)
colnames(run2_complete) <- c("condition",rep("WT",30))
head(run2_complete)

# write csv
write.csv(run2_complete, file="two_batches_run2_for_CSI.csv",row.names = F)



###### Extract data for three batch approach: ######

three_batch_runA <- read.csv("three_batch_runA.csv", header=F)
head(three_batch_runA)
three_batch_runB <- read.csv("three_batch_runB.csv", header=F)
head(three_batch_runB)
three_batch_runC <- read.csv("three_batch_runC.csv", header=F)
head(three_batch_runC)

# select only genes which are in the TF families selected for three_batch_runA
runA <- tpm_data[which(tpm_data$family %in% three_batch_runA$V1),1:31]
head(runA)
dim(runA)
#add correct column names
colnames(runA) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(runA)
runA <- rbind(colnames(runA), runA)
colnames(runA) <- c("condition",rep("WT",30))

head(runA)
write.csv(runA, file="three_batches_runA_for_CSI.csv",row.names = F)


# select only genes which are in the TF families selected for three_batch_runB
runB <- tpm_data[which(tpm_data$family %in% three_batch_runB$V1),1:31]
head(runB)
dim(runB)
#add correct column names
colnames(runB) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(runB)
runB <- rbind(colnames(runB), runB)
colnames(runB) <- c("condition",rep("WT",30))

head(runB)
write.csv(runB, file="three_batches_runB_for_CSI.csv",row.names = F)

# select only genes which are in the TF families selected for three_batch_runC
runC <- tpm_data[which(tpm_data$family %in% three_batch_runC$V1),1:31]
head(runC)
dim(runC)
#add correct column names
colnames(runC) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(runC)
runC <- rbind(colnames(runC), runC)
colnames(runC) <- c("condition",rep("WT",30))

head(runC)
write.csv(runC, file="three_batches_runC_for_CSI.csv",row.names = F)
