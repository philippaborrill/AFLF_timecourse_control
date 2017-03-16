#!/bin/bash
#
# SLURM batch script to launch BLAST
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/extractNACbyBLAST.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/extractNACbyBLAST.%N.%j.err # STDERR
#SBATCH -J extractNAC
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/


source blast+-2.2.28

#makeblastdb -in Tritcum_aestivum_seq.fas -parse_seqids -dbtype prot

#grep "|NAC|" Tritcum_aestivum_seq.fas > Triticum_aestivum_NAC_IDs.txt

#blastdbcmd -db Tritcum_aestivum_seq.fas -dbtype prot -entry_batch Triticum_aestivum_NAC_IDs.txt -out NAC_Triticum_aestivum.fa

blastdbcmd -db Tritcum_aestivum_seq.fas -dbtype prot -entry Tae000413|Triticum_aestivum|NAC|PUT-163b-Triticum_aestivum-46870;gnl|UG|Ta#S61799989;PUT-163b-Triticum_aestivum Triticum_aestivum_NAC_IDs.txt -out NAC_Triticum_aestivum.fa -46869;gnl|UG|Ta#S24917174;gnl|UG|Ta#S61779460
