# script to select genes for CSI for 10 tpm in specific families

# Philippa Borrill
# 30.09.2016


setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\9_EBSeqHMM_TF_only_new_annotation\\")

# read file which was produced by the scipt TGAC_EBSeqHMM_TF_only_try_different_tpm_and_prep_data_TCAP_CSI.R
tpm_data <- read.csv("countData_diff_expr_TF_10tpm_high_conf_with_TF.csv", header=T)
head(tpm_data)

# now read in file containing list of NAC, MYB and WRKY which were found in CSI network 10 tpm 0.045 threshold with new CSI running
TF_in_network <- read.csv("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\8_CSI_analysis\\new_annotation_NAC_WRKY_MYB\\new_CSI_new_annotation\\NAC_MYB_WRKY_10tpm\\0.045_network_small_v.csv", header=T)
head(TF_in_network)




###### Extract data for two batch approach: ######

setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\8_CSI_analysis\\new_annotation_NAC_WRKY_MYB\\new_CSI_new_annotation\\other_families")


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





###### Extract data for four batch approach: ######

setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\8_CSI_analysis\\new_annotation_NAC_WRKY_MYB\\new_CSI_new_annotation\\other_families")

four_batch_run1 <- read.csv("four_batch_run1.csv", header=F)
head(four_batch_run1)
four_batch_run2 <- read.csv("four_batch_run2.csv", header=F)
head(four_batch_run2)
four_batch_run3 <- read.csv("four_batch_run3.csv", header=F)
head(four_batch_run3)
four_batch_run4 <- read.csv("four_batch_run4.csv", header=F)
head(four_batch_run4)

# select only genes which are in the TF_network from CSI
data_TF_in_network <- merge(tpm_data, TF_in_network, by.x = "condition", by.y = "name")
head(data_TF_in_network)
colnames(data_TF_in_network)[1:31]
data_TF_in_network <- data_TF_in_network[,1:31]
head(data_TF_in_network)
dim(data_TF_in_network)

# select only genes which are in the TF families selected for four_batch_run1
colnames(tpm_data)

run1 <- tpm_data[which(tpm_data$family %in% four_batch_run1$V1),1:31]
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
write.csv(run1_complete, file="four_batches_run1_for_CSI.csv",row.names = F)


# select only genes which are in the TF families selected for four_batch_run2
run2 <- tpm_data[which(tpm_data$family %in% four_batch_run2$V1),1:31]
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
write.csv(run2_complete, file="four_batches_run2_for_CSI.csv",row.names = F)

# select only genes which are in the TF families selected for four_batch_run3
run3 <- tpm_data[which(tpm_data$family %in% four_batch_run3$V1),1:31]
head(run3)
dim(run3)

# add together the NAC,MYB,WKRY data with the run3 data
run3_complete <- rbind.data.frame(data_TF_in_network, run3)
head(run3_complete)
dim(run3_complete)

# add in correct column names
colnames(run3_complete) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(run3_complete)
run3_complete <- rbind(colnames(run3_complete), run3_complete)
colnames(run3_complete) <- c("condition",rep("WT",30))
head(run3_complete)

# write csv
write.csv(run3_complete, file="four_batches_run3_for_CSI.csv",row.names = F)


# select only genes which are in the TF families selected for four_batch_run4
run4 <- tpm_data[which(tpm_data$family %in% four_batch_run4$V1),1:31]
head(run4)
dim(run4)

# add together the NAC,MYB,WKRY data with the run4 data
run4_complete <- rbind.data.frame(data_TF_in_network, run4)
head(run4_complete)
dim(run4_complete)

# add in correct column names
colnames(run4_complete) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(run4_complete)
run4_complete <- rbind(colnames(run4_complete), run4_complete)
colnames(run4_complete) <- c("condition",rep("WT",30))
head(run4_complete)

# write csv
write.csv(run4_complete, file="four_batches_run4_for_CSI.csv",row.names = F)


###### now make supernetwork input files for CSI #####
setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\8_CSI_analysis\\new_annotation_NAC_WRKY_MYB\\all_families")

supernetwork1 <- read.csv("supernetwork1_all_genes_run1_to_run4_threshold_0.013.csv", header=T)
head(supernetwork1)
dim(supernetwork1)

# select only genes which are in the list from CSI (all genes run 1 to run 4)
data_TF_in_network <- merge(tpm_data, supernetwork1, by.x = "condition", by.y = "genes")
head(data_TF_in_network)
colnames(data_TF_in_network)[1:31]
data_TF_in_network <- data_TF_in_network[,1:31]
head(data_TF_in_network)
dim(data_TF_in_network)

# add in correct column names
colnames(data_TF_in_network) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(data_TF_in_network)
data_TF_in_network <- rbind(colnames(data_TF_in_network), data_TF_in_network)
colnames(data_TF_in_network) <- c("condition",rep("WT",30))
head(data_TF_in_network)

# write csv
write.csv(data_TF_in_network, file="supernetwork1_for_CSI.csv",row.names = F)


## do it for supernetwork2

supernetwork2 <- read.csv("supernetwork2_all_genes_run1_to_run4_threshold_0.013_2_degrees.csv", header=T)
head(supernetwork2)
dim(supernetwork2)

# select only genes which are in the list from CSI (all genes run 1 to run 4)
data_TF_in_network <- merge(tpm_data, supernetwork2, by.x = "condition", by.y = "genes")
head(data_TF_in_network)
colnames(data_TF_in_network)[1:31]
data_TF_in_network <- data_TF_in_network[,1:31]
head(data_TF_in_network)
dim(data_TF_in_network)

# add in correct column names
colnames(data_TF_in_network) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(data_TF_in_network)
data_TF_in_network <- rbind(colnames(data_TF_in_network), data_TF_in_network)
colnames(data_TF_in_network) <- c("condition",rep("WT",30))
head(data_TF_in_network)

# write csv
write.csv(data_TF_in_network, file="supernetwork2_for_CSI.csv",row.names = F)


##### extract data to run supernetwork of all TF found in 4 batches at threshold 0.05 13/10/2016 ##### 

setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\9_EBSeqHMM_TF_only_new_annotation\\")

# read file which was produced by the scipt TGAC_EBSeqHMM_TF_only_try_different_tpm_and_prep_data_TCAP_CSI.R
tpm_data <- read.csv("countData_diff_expr_TF_10tpm_high_conf_with_TF.csv", header=T)
head(tpm_data)

setwd("Y:\\PB_AFLF\\control_timecourse\\TGAC_kallisto_analysis\\kallisto_results_bootstrap\\results\\8_CSI_analysis\\new_annotation_NAC_WRKY_MYB\\new_CSI_new_annotation\\other_families\\four_batch")

supernetwork <- read.csv("list_of_all_genes_from_4_runs.csv", header=F)
head(supernetwork)
dim(supernetwork)

# select only genes which are in the list from CSI (all genes run 1 to run 4)
data_TF_in_network <- merge(tpm_data, supernetwork, by.x = "condition", by.y = "V1")
head(data_TF_in_network)
colnames(data_TF_in_network)[1:31]
data_TF_in_network <- data_TF_in_network[,1:31]
head(data_TF_in_network)
dim(data_TF_in_network)

# add in correct column names
colnames(data_TF_in_network) <- c("time",rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),each=3))
head(data_TF_in_network)
data_TF_in_network <- rbind(colnames(data_TF_in_network), data_TF_in_network)
colnames(data_TF_in_network) <- c("condition",rep("WT",30))
head(data_TF_in_network)

# write csv
write.csv(data_TF_in_network, file="supernetwork_for_CSI_4_runs_new_CSI_new_annot_10tpm.csv",row.names = F)

