#!/usr/bin/perl -w

# 20-03-2017
# Philippa.borrill@jic.ac.uk
#
# Aim of script is to add back in full gene IDs to tree (Raxml only works with 10 character phylip files)

use File::Copy;
my $path = '/nbi/Research-Groups/ifs/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TF_analysis/phylogenetics';
my $output_tree = 'RAxML_bestTree.Triticum_aestivum_NACs_msa_tree_copy_longIDs';

chdir("$path") or die "couldn't move to directory";

my %dict;
open(GENEIDS, "geneIDs.txt");
while(<GENEIDS>){
	/(\S+)\s+(\S+)/;
$dict{$1}={$2};
}
close GENEIDS;

open(TREEFILE, "RAxML_bestTree.Triticum_aestivum_NACs_msa_tree_copy");
while(my $line = <>){
    foreach my $s (keys %dict){
        $line =~ s/$s/$dict{$s}/g;
    }
		open (OUTPUT_FILE, ">>$output_tree") or die "Couldn't open $output_tree file\n";
		print OUTPUT_FILE "$line";
		close OUTPUT_FILE;
}

close TREEFILE;
