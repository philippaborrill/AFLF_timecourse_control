#!/bin/bash
#
# SLURM script to launch WGCNA in R
#
#SBATCH -p nbi-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem 250000 # memory pool for all cores (MB)
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/slurm_output/WGCNA.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/slurm_output/WGCNA.%N.%j.err # STDERR
#SBATCH -J WGCNA
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address


cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/scripts/WGCNA

source R-3.1.0

Rscript part2_WGCNA_expVIP_cluster.R
