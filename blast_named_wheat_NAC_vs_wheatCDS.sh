#!/bin/bash
#
# SLURM batch script to launch BLAST
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/BLAST_wheat_NAC_wheat.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/BLAST_wheat_NAC_wheat.%N.%j.err # STDERR
#SBATCH -J blast
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

#cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/table_orthologues
#source blast+-2.2.28
#makeblastdb -in Triticum_aestivum.TGACv1.cds.all.fa -dbtype nucl -parse_seqids

source blast+-2.2.28
cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/table_orthologues
blastn -db Triticum_aestivum.TGACv1.cds.all.fa -query wheat_NACs_with_names.txt -num_threads 1  -max_target_seqs 1 -outfmt 6 -out wheat_NACs_with_names_blast_to_wheat.txt
