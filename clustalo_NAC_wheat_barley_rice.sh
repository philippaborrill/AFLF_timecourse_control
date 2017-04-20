#!/bin/bash
#
# SLURM batch script to clustalo_NAC_Triticum_aestivum
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/clustalo_NAC_wheat_barley_rice.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/clustalo_NAC_wheat_barley_rice.%N.%j.err # STDERR
#SBATCH -J clustalo_NAC_wheat_barley_rice
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/phylogenetics

source clustalo-1.2.0

#clustalo -i Triticum_aestivum_NACs_longest_isoforms.fa -o Triticum_aestivum_NACs_msa.clu2 --outfmt=clu
clustalo -i wheat_barley_rice_NAC_longest_isoforms.fa -o wheat_barley_rice_NAC_msa.clu --outfmt=clu

#if want fasta format
clustalo -i wheat_barley_rice_NAC_longest_isoforms.fa -o wheat_barley_rice_NAC_msa.fa
