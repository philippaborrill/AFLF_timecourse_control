#!/bin/bash
#
# SLURM batch script to launch BLAST
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/BLAST.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/BLAST.%N.%j.err # STDERR
#SBATCH -J blast
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

#cd /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release
#source blast+-2.2.28
#makeblastdb -in Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.pep_no_tabs.fa -dbtype prot -parse_seqids

source blast+-2.2.28
cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis
#blastp -db /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.pep_no_tabs.fa -query Tritcum_aestivum_seq.fas -num_threads 1  -max_target_seqs 10 -outfmt 6 -out Tritcum_aestivum_blast_to_TGAC.txt

#blastp -db /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.pep_no_tabs.fa -query Aegilops_tauschii_seq.fas -num_threads 1  -max_target_seqs 10 -outfmt 6 -out Aegilops_tauschii_blast_to_TGAC.txt

#blastp -db /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.pep_no_tabs.fa -query Hordeum_vulgare_seq.fas -num_threads 1  -max_target_seqs 10 -outfmt 6 -out Hordeum_vulgare_blast_to_TGAC.txt

#blastp -db /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.pep_no_tabs.fa -query Oryza_indica_seq.fas -num_threads 1  -max_target_seqs 10 -outfmt 6 -out Oryza_indica_blast_to_TGAC.txt

#blastp -db /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.pep_no_tabs.fa -query Oryza_japonica_seq.fas -num_threads 1  -max_target_seqs 10 -outfmt 6 -out Oryza_japonica_blast_to_TGAC.txt

blastp -db /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.pep_no_tabs.fa -query Triticum_urartu_seq.fas -num_threads 1  -max_target_seqs 10 -outfmt 6 -out Triticum_urartu_blast_to_TGAC.txt
