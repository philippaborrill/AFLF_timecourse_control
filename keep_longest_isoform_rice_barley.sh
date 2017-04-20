#!/bin/bash
#
# SLURM batch script to extract longest isoforms
# from https://www.biostars.org/p/107759/
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/extract_longest-isoform_rice_barley.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/slurm_output/extract_longest-isoform_rice_barley.%N.%j.err # STDERR
#SBATCH -J extract_long_isoform
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/phylogenetics

cat hordeum_vulgare_NAC.fasta | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | tr "." "\t" | awk -F '	'  '{printf("%s\t%d\n",$0,length($3));}' | sort -t '	' -k1,1 -k4,4nr | sort -t '	' -k1,1 -u -s | sed 's/	/./' | cut -f 1,2 | tr "\t" "\n"  | fold -w 80 > hordeum_vulgare_NAC_longest_isoforms.fa
cat oryza_sativa_japonica_NAC.fasta | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | tr "." "\t" | awk -F '	'  '{printf("%s\t%d\n",$0,length($3));}' | sort -t '	' -k1,1 -k4,4nr | sort -t '	' -k1,1 -u -s | sed 's/	/./' | cut -f 1,2 | tr "\t" "\n"  | fold -w 80 > oryza_sativa_japonica_NAC_longest_isoforms.fa


#COMMENTS:
#cat Triticum_aestivum_NACs.fa |\
#awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' |\ #linearize fasta
#tr "." "\t" |\ #extract version from header
#awk -F '	'  '{printf("%s\t%d\n",$0,length($3));}' |\ #extact length
#sort -t '	' -k1,1 -k4,4nr |\ #sort on name, inverse length
#sort -t '	' -k1,1 -u -s |\ #sort on name, unique, stable sort (keep previous order)
#sed 's/	/./' |\ #restore name
#cut -f 1,2 |\ #cut name, sequence
#tr "\t" "\n"  |\ #back to fasta
#fold -w 60 |\ #pretty fasta
#> Triticum_aestivum_NACs_longest_isoforms.fa #write output file
