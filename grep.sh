#!/bin/bash
#
# SLURM batch script to launch grep
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/BLAST.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/BLAST.%N.%j.err # STDERR
#SBATCH -J grep
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis
grep -f unique_TFs_BLAST_ensembl_GO.txt /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.pep.fa >> unique_TFs_BLAST_ensembl_GO_full_ID.txt

