#!/bin/bash
#
# SLURM batch script to launch BLAST
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/BLAST_wheat_NAC_rice.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/BLAST_wheat_NAC_rice.%N.%j.err # STDERR
#SBATCH -J blast
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

#cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/table_orthologues
#source blast+-2.2.28
#makeblastdb -in o_sativa_v7.0_MSU_cds.fa -dbtype nucl -parse_seqids

source blast+-2.2.28
cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/table_orthologues
blastn -db o_sativa_v7.0_MSU_cds.fa -query wheat_high_conf_NAC_CDS.txt -num_threads 1  -max_target_seqs 1 -outfmt 6 -out wheat_high_conf_NAC_CDS_blast_to_rice.txt
