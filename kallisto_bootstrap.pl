#!/usr/bin/perl -w

# Re-engineered to use LSF properly
# Philippa.borrill@jic.ac.uk
#
# Aim of script is to run kallisto on RNA-seq for multiple samples to a common reference to calculate expression levels

#### paths and references:
my $path = '/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/';
my $read_path = "PB_AFLF/control_timecourse/raw_data/PKG_ENQ-789_60_samples_data_transfer/";
my $ref = "/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/IWGSC_data/Triticum_aestivum.IWGSC2.26.cdna.all.fa";
my $index = "/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/expression_browser/kallisto_analysis/Triticum_aestivum.IWGSC2.26.cdna.all";

# NB make index by kallisto index -i Triticum_aestivum.IWGSC2.26.cdna.all Triticum_aestivum.IWGSC2.26.cdna.all.fa

my $output_dir = "$path/PB_AFLF/control_timecourse/kallisto_analysis/";

### input info: contains 5 tab separated columns with: directory, R1, R2, R1_pair2, R2_pair2, sample_name
### must be in $output_dir
my $control_RNA_seq_paired_list = 'input_for_kallisto.txt';

#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$output_dir") or die "couldn't move to output directory";

open (INPUT_FILE, "$control_RNA_seq_paired_list") || die "couldn't open the input file $control_RNA_seq_paired_list!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";
			
#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";

my $dir = $array[0];
my $pair_1_R1 = $array[1];
my $pair_1_R2 = $array[2];
my $pair_2_R1 = $array[3];
my $pair_2_R2 = $array[4];
my $output = $array[5];

#print "$path/$read_path/$dir\n";

chdir("$path/$read_path/$dir") or die "couldn't move to specific read directory";

### bsub header including memory usage request
my $bsub_header = <<"LSF";
#!/bin/bash
#
# LSF batch script to launch parallel kallisto tasks
#
#BSUB -q NBI-Prod128
#BSUB -J kallisto_JOBNAME
#BSUB -R "rusage[mem=10000]"
#BSUB -n 8
LSF


 my $tmp_file = "$output_dir/tmp/samtools_kallisto_paired_lsf.$dir";

  open (BSUB, ">$tmp_file") or die "Couldn't open temp file\n";
  $bsub_header =~ s/JOBNAME/$output/;
  print BSUB "$bsub_header\n\n";
  print BSUB "\ncd $path/$read_path/$dir\n";

  print BSUB "source kallisto-0.42.3\n";
  print BSUB "source samtools-0.1.19\n";

#	print BSUB "kallisto quant -i $index -o $output_dir/$output -b 30 -t 8 --pseudobam $pair_1_R1 $pair_1_R2 $pair_2_R1 $pair_2_R2 | samtools view -Sb - > $output_dir/$output".".bam \n";
  print BSUB "samtools sort $output_dir/$output".".bam $output_dir/$output".".sorted\n";
  print BSUB "samtools index $output_dir/$output".".sorted.bam $output_dir/$output".".sorted.bai\n";

  close BSUB;
  system("bsub < $tmp_file");
 # unlink $tmp_file;

## need to close loop which goes through all of the directories in the list
	}
	    close(INPUT_FILE); 

#BSUB -o $output_dir/kallisto_bootstrap_log.txt

