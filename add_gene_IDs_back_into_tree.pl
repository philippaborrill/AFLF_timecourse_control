#!/usr/bin/perl -w

# 20-03-2017
# Philippa.borrill@jic.ac.uk
#
# Aim of script is to add back in full gene IDs to tree (Raxml only works with 10 character phylip files)

use File::Copy;
my $path = '/nbi/Research-Groups/ifs/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/phylogenetics';

my $ID_file = 'geneIDs.txt';
# $ID_file must have two columns (tab separated) - first with abbreviated gene IDs
# e.g. AA1047550 and 2nd with full gene ID e.g. TRIAE_CS42_4BL_TGACv1_320737_AA1047550


my $tree = 'RAxML_bestTree.Triticum_aestivum_NACs_msa_tree_copy';

my $output_tree = 'RAxML_bestTree.Triticum_aestivum_NACs_msa_tree_copy_longIDs';

chdir("$path") or die "couldn't move to directory";

# open ID_file to get shortID and longID which will be used for substitution
		open (INPUT_FILE2, "$ID_file") || die "couldn't open the input file $ID_file!";
				    while (my $line2 = <INPUT_FILE2>) {
					chomp $line2;
					my @array = split(/\t/,$line2);

					my $shortID = $array[0];
					my $longID = $array[1];

# open $tree file which needs to be edited
					open (INPUT_FILE1, "$tree") || die "couldn't open the input file $tree!";
							    while (my $line = <INPUT_FILE1>) {
										chomp $line;
# substitute $shortID for $longID in the $line
										$line =~ s/$shortID/$longID/g;

# print the edited $line to $output_tree file
										open (OUTPUT_FILE, ">$output_tree") or die "Couldn't open $output_tree file\n";
										print OUTPUT_FILE "$line";
										close OUTPUT_FILE;
										# copy $output_tree (edited) to $tree to re-use as input
										copy($output_tree, $tree) or die "copy failed: $!";
											}
			close(INPUT_FILE1);

# need to close loop which goes through all of the gene IDs in the list
	}
	    close(INPUT_FILE2);
