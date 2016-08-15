#!/bin/bash
#
# SLURM script to launch kallisto index
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/slurm_output/kallisto_index.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/slurm_output/kallisto_index.%N.%j.err # STDERR
#SBATCH -J kallisto_index
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address

cd  /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release

source kallisto-0.42.3

kallisto index -i Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.cdna Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.cdna.fa
