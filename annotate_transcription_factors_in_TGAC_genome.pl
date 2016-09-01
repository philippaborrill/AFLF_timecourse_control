#!/nbi/software/production/perl/5.16.2/x86_64/bin/perl/ -w

use strict;
use diagnostics;
use warnings;

# Aim of this script is to take a file with TGAC gene names and GO and PFAM annotations and find genes which have specific PFAM annotations (related to being transcription factors). These transcription factor (TF) PFAMs are stored in a separate file. A new file is create with the original file's information + 1 new column saying the PFAM if it is a transcription factor + 1 column with the family of the transcription factor

my $path = "/nbi/group-data/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/1_SOM_analysis";

#my $genome_file = "/nbi/group-data/NBI/Cristobal-Uauy/PB_AFLF/control_timecourse/TGAC_kallisto_analysis/kallisto_results_bootstrap/results/1_SOM_analysis/test_annotation.tsv";
my $genome_file = "/nbi/group-data/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.functional_annotation.tsv"; # file with all TGAC genes with PFAM, GO term information, has 17 columns (tab separated)

# to run programme need to provide a text file containing a list of transcription factor families with PFAM terms associated - must be saved in $path
#my $TF_PFAM_file = "test_transcription_factor_family_PFAM.txt";
my $TF_PFAM_file = "transcription_factor_family_PFAM.txt";


my $outputfile = "TGAC_genes_with_$TF_PFAM_file";


#change to directory containing the text files
chdir ("$path") or die;


open (OUTPUTFILE, ">>$outputfile") or die "Cannot open file $outputfile to write to\n\n";
print OUTPUTFILE "Transcription factor PFAM file was: $TF_PFAM_file\n";
print OUTPUTFILE "Gene list file was: $genome_file\n\n ";
print OUTPUTFILE "#1.Transcript_ID\t#2.Gene_ID\t#3.Gene_biotype\t#4.Representative\t#5.Gene_confidence\t#6.Transcript_confidence\t#7.Transposon_related\t#8.AHRD_quality_code\t#9.AHRD_description\t#10.AHRD_GO_term\t#11.PFAM_ID\t#12.Interproscan_family\t#13.Interproscan_domain\t#14.Interproscan_sites/repeats\t#15.Interproscan_GO_term\t#16.Interproscan_EC_number\t#17.Interproscan_pathways\t#18.TF_PFAM\t#19.TF_family\n ";
close (OUTPUTFILE);

#open the text file
open(TF,"$TF_PFAM_file") or die "Could not open $TF_PFAM_file: $!";

while (my $line = <TF>){
		chomp $line;
    	my @array = split(/\t/,$line);

my $family = $array[0];
my $PFAM = $array[1];
print "$PFAM\n";


	open(GENE_LIST, "$genome_file") or die "Could not open $genome_file: $!";
	while (my $gene_line = <GENE_LIST>){
print "$gene_line\n";
chomp $gene_line;
		
# if PFAM doesn't match line in genome_file print genome_file original line
# if PFAM does match line in genome_file print genome_file original line and PFAM and family name

		if ($gene_line =~ /$PFAM/ ) { 
				open(OUTPUT, ">>$outputfile") or die "Cannot open file $outputfile to write to\n\n";
			    	print OUTPUT "$gene_line\t$PFAM\t$family\n";
				close(OUTPUT);	
				}		
	}
}
