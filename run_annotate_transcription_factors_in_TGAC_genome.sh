#!/bin/bash
#
# SLURM batch script to launch annotate_transcription_factors_in_TGAC_genome.pl
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 10000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/group-data/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/slurm_output/annotate_TFs.%N.%j.out # STDOUT
#SBATCH -e /nbi/group-data/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/slurm_output/annotate_TFs.%N.%j.err # STDERR
#SBATCH -J annotate_TFs
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address


cd /nbi/group-data/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/scripts

perl annotate_transcription_factors_in_TGAC_genome.pl
