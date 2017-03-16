#!/bin/bash
#
# SLURM batch script to raxml_NAC_Triticum_aestivum
# used http://evomics.org/learning/phylogenetics/raxml/ as a guide
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/raxml_NAC_Triticum_aestivum.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/raxml_NAC_Triticum_aestivum.%N.%j.err # STDERR
#SBATCH -J raxml_NAC_Triticum_aestivum
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/phylogenetics

source raxml-8.1.2.x


raxmlHPC -s Triticum_aestivum_NACs_msa.fa -n Triticum_aestivum_NACs_msa_tree -m PROTGAMMAAUTO -f a -x 100 -N 100 -p 100
