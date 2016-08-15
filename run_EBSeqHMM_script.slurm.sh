#!/bin/bash
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 250000 # memory pool for all cores (MB) # this job needs lots of memory
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/scripts/EBSeqHMM10_02_2016.%N.%j.out # STDOUT
#SBATCH -e /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/scripts/EBSeqHMM10_02_2016.%N.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address



cd /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/scripts/

source R-3.1.0

Rscript EBSeqHMM_script.R

