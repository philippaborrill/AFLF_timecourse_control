#!/usr/bin/perl -w

# Philippa.borrill@jic.ac.uk
# 08-05-17
#
# Aim of script is to run RaXML on different NAC subgroups

#### paths and references:
my $path = '/nbi/Research-Groups/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/phylogenetics/trees/groups_a-h/';

## input info is a list of files to run RaXML on (PHYLIP 4 format)
my $list = 'list_of_group_files.txt';

#############################

#open the input file and go through the lines one by one
chdir("$path") or die "couldn't move to directory";

open (INPUT_FILE, "$list") || die "couldn't open the input file $list!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;

			my $tmp_file = "$path/tmp/RaxML.$line";

my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel RaxML tasks
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 3000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o $path/slurm_output/raxml.JOBNAME.%N.%j.out # STDOUT
#SBATCH -e $path/slurm_output/raxml.JOBNAME.%N.%j.err # STDERR
#SBATCH -J JOBNAME_RaxML
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address
SLURM

		  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
		   $SLURM_header =~ s/JOBNAME/$line/g;

			 print SLURM "$SLURM_header\n\n";
			  print SLURM "cd $path\n";
			  print SLURM "source raxml-8.1.2.x\n";
   			print SLURM "raxmlHPC -s $line -n $line.tree -m PROTGAMMAAUTO -f a -x 100 -N 100 -p 100 \n";

			  close SLURM;
			  system("sbatch $tmp_file");
			 # unlink $tmp_file;

			## need to close loop which goes through all of the directories in the list
				}
				    close(INPUT_FILE);
