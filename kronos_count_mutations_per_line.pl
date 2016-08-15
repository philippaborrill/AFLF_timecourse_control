#!/nbi/software/production/perl/5.16.2/x86_64/bin/perl/ -w

use strict;
use diagnostics;
use warnings;

# Aim of this script is to take a file with gene names in a column and find all mutations in those genes. It makes a new file containing only these genes/mutations. It then counts the number of mutations per gene, the number of stop/splice/start codon mutations per gene and outputs this data to a table. 

my $path = "/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/TILLING/";

my $dir = "/count_mutations/";


my $mutation_file = "$path/$dir/data/Kronos_June_2015-SIFT-final.tsv"; # kronos mutation file 

# to run programme need to provide a text file containing a list of directories
unless ($ARGV[0]) {
  print "Usage: kronos_count_mutants.pl <text file> (text file containing list of genes - must saved be in $path/$dir)\n";
  exit;
}


my $input_file = $ARGV[0];

my $outputfile = "kronos_list_of_mutations_$input_file";

my $count_file = "kronos_count_of_mutations_$input_file";

#change to directory containing the text files
chdir ("$path/$dir") or die;

#first print to final output files what the inputs were:
open (COUNTFILE, ">>$count_file") or die "Cannot open file $count_file to write to\n\n";
print COUNTFILE "Inputfile was: $input_file\n";
print COUNTFILE "Mutation file was: $mutation_file\n ";
print COUNTFILE "Line\tTotal_mutations\tTruncations\n"; 
close (COUNTFILE);

open (OUTPUTFILE, ">>$outputfile") or die "Cannot open file $outputfile to write to\n\n";
print OUTPUTFILE "Inputfile was: $input_file\n";
print OUTPUTFILE "Mutation files was: $mutation_file\n\n ";
print OUTPUTFILE "Scaffold\tChromosome\tLibrary\tLine\tPosition\tChromosome_Start\tRef_base\tWT_base\tmut_base\thet/hom\tWT_cov\tmut_cov\tconfidence\tGene_model\tTranscript_model\tConsequence\tcDNA_pos\tCDS_pos\tAmino_acids\tCodons\tSIFT\n ";
close (OUTPUTFILE);

#open the text file
open(GENES,"$input_file") or die "Could not open $input_file: $!";

while (my $line = <GENES>){
	print $line;
	chomp $line;
    	
	open(MUTATIONS, "$mutation_file") or die "Could not open $mutation_file: $!";
	while (my $mutation = <MUTATIONS>){
		
# print each mutation to output file
		if ($mutation =~ /Kronos$line/) {
				open(OUTPUT, ">>$outputfile") or die "Cannot open file $outputfile to write to\n\n";
			    	print OUTPUT "$mutation";
				close(OUTPUT);	
			}		
		}


### steps to count total number and number of truncation mutations
# set parameters
my $total_mutations= "0";
my $stop_mutations= "0";
my $splice_donor_mutations= "0";
my $splice_acceptor_mutations= "0";
my $start_lost_mutations= "0";
my $initiator_codon_mutations= "0";

open(OUTPUT, '<', "$outputfile") or die "Cannot open file $outputfile\n";

while (my $row = <OUTPUT>){
		if ($row =~ /$line/) {
				$total_mutations = $total_mutations+1;
				}

		if ($row =~ /$line/ and $row =~ /stop_gained/) {
				$stop_mutations = $stop_mutations+1;
				}

		if ($row =~ /$line/ and $row =~ /splice_donor_variant/) {
				$splice_donor_mutations = $splice_donor_mutations+1;
				}

		if ($row =~ /$line/ and $row =~ /splice_acceptor_variant/) {
				$splice_acceptor_mutations = $splice_acceptor_mutations+1;
				}

		if ($row =~ /$line/ and $row =~ /start_lost/) {
				$start_lost_mutations = $start_lost_mutations+1;
				}

		if ($row =~ /$line/ and $row =~ /initiator_codon_variant/) {
				$initiator_codon_mutations = $initiator_codon_mutations+1;
				}
		}

close(OUTPUT);
	
# step to add up all truncation mutations in each gene
	my $truncation_mutations = $stop_mutations + $splice_donor_mutations + $splice_acceptor_mutations + $start_lost_mutations + $initiator_codon_mutations;

	# write the number of total and trucation mutations for that gene to file
	open(COUNTFILE, ">>$count_file") or die "Cannot open file $count_file to write to\n\n";
	print COUNTFILE "$line\t$total_mutations\t$truncation_mutations\n";
	close (COUNTFILE);
		
}
	
