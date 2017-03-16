#!/bin/bash
#
# SLURM batch script to launch BLAST
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/extractNAC.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/extractNAC.%N.%j.err # STDERR
#SBATCH -J extractNAC
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

source perl-5.22.1

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/
perl /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/scripts/fasta_extract_by_pattern.pl -r -p "|NAC|" Tritcum_aestivum_seq.fas > Triticum_aestivum_NAC.fa