#!/bin/bash
#
# SLURM batch script to merge results
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/kallisto_results_bootstrap/slurm_output/merge_samples.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/kallisto_results_bootstrap/slurm_output/merge_samples.%N.%j.err # STDERR
#SBATCH -J merge_samples
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address

##############
# usage: make sure you have the previous datasets in $MY_PREVIOUS_RESULTS_DIRECTORY (keep a back-up copy somewhere else too!)
# these should be files: edited_final_output_counts.txt edited_final_output_tpm.txt
# make sure new samples to add (folders containing eXpress results) are in directory: $MY_WORKING_DIRECTORY

######################## functions #################

# body will allow you to do a function on a text file excluding the header row
# print the header (the first line of input)
# and then run the specified command on the body (the rest of the input)
# use it in a pipeline, e.g. ps | body grep somepattern
function body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}
####################################################

##### move the NEW results dir e.g. SRR1958738 to a new dir which must be $MY_WORKING_DIRECTORY
MY_WORKING_DIRECTORY=/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly_Jemima/analysis/kallisto_results_bootstrap/results/
MY_PREVIOUS_RESULTS_DIRECTORY=/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/kallisto_results_bootstrap/results/

# move to top directory  
cd $MY_WORKING_DIRECTORY

# need file with list of gene names in 1 column - in sorted order - as defined by unix with header: transcript
MY_TRANSCRIPTS_ORDERED=/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/kallisto_results_bootstrap/results/transcripts_ordered.txt

# make empty output + temp files:
echo "" > $MY_WORKING_DIRECTORY/output_counts.txt
echo "" > $MY_WORKING_DIRECTORY/output_tpm.txt
echo "" > $MY_WORKING_DIRECTORY/temp.txt
echo "" > $MY_WORKING_DIRECTORY/col4.txt
echo "" > $MY_WORKING_DIRECTORY/col5.txt

#  two output files: output_counts.txt and output_tpm.txt
MY_OUTPUT_COUNTS=$MY_WORKING_DIRECTORY/output_counts.txt

MY_OUTPUT_TPM=$MY_WORKING_DIRECTORY/output_tpm.txt

#temp file:
MY_TEMP_FILE=$MY_WORKING_DIRECTORY/temp.txt

#col4 and 5 files:
MY_COL4=$MY_WORKING_DIRECTORY/col4.txt
MY_COL5=$MY_WORKING_DIRECTORY/col5.txt

# put the list of ordered gene names as the left hand column in the output files:
cp $MY_TRANSCRIPTS_ORDERED $MY_OUTPUT_COUNTS
cp $MY_TRANSCRIPTS_ORDERED $MY_OUTPUT_TPM

# loop through each dir ($dir_name) in turn which has a results file in it and put the values of counts and tpm into relevant output files

	for dir_name in */ ; do
	#for dir_name in `ls -d */`; do
	
   	echo "directory name is: $dir_name \n"
	
	
	cd $dir_name
	# check abundance.tsv has 273,740 lines
	MY_WC=`cat abundance.tsv | wc -l`
	echo "My line count is: $MY_WC \n"
	
		
		if [ "$MY_WC" == "273740" ]

		then
		# sort abundance.tsv on the 1st column

		cat abundance.tsv | body sort -k1,1 > sorted.abundance.tsv


		# replace headings of 4th and 5th column with name of directory

		sed "s|est_counts|$dir_name|" sorted.abundance.tsv > renamed.sorted.abundance.tsv
		sed "s|tpm|$dir_name|" renamed.sorted.abundance.tsv > renamed2.sorted.abundance.tsv
		mv renamed2.sorted.abundance.tsv renamed.sorted.abundance.tsv
		
		# put 4th column into counts output, 5th column into tpm output
		# 1st need to cut then paste, then rename the paste file back to be the new input file:
		cut -f 4 renamed.sorted.abundance.tsv > $MY_COL4	
		cut -f 5 renamed.sorted.abundance.tsv > $MY_COL5

		#empty $MY_TEMP_FILE
		echo "" > $MY_TEMP_FILE		
		paste $MY_OUTPUT_COUNTS $MY_COL4 > $MY_TEMP_FILE
		mv $MY_TEMP_FILE $MY_OUTPUT_COUNTS
		
			
		#empty $MY_TEMP_FILE
		echo "" > $MY_TEMP_FILE	
		paste $MY_OUTPUT_TPM $MY_COL5 > $MY_TEMP_FILE
		mv $MY_TEMP_FILE $MY_OUTPUT_TPM
	
		else 
		echo "The abundance.tsv file in $dir_name did not have the correct number of lines therefore was ignored"
		fi

		cd ..
	#next dir
# end loop
	done

# replace / at end of each directory
sed "s|/||g" $MY_OUTPUT_COUNTS > final_output_counts.txt
sed "s|/||g" $MY_OUTPUT_TPM > final_output_tpm.txt

# add step to merge new samples with previous samples

######## don't need this currently so write exit

# previous samples are stored in: $MY_PREVIOUS_RESULTS_DIRECTORY (see at the top)

# have 2 files there:
MY_PREVIOUS_OUTPUT_COUNTS=$MY_PREVIOUS_RESULTS_DIRECTORY/final_output_counts.txt
MY_PREVIOUS_OUTPUT_TPM=$MY_PREVIOUS_RESULTS_DIRECTORY/final_output_tpm.txt

# remove first column from new data sheets
cut -f 2- final_output_counts.txt > cut_final_output_counts.txt
cut -f 2- final_output_tpm.txt > cut_final_output_tpm.txt

# combine the new and old data sheets in the folder "working_dir" and renamed as edited_final_output...txt (so this script can be re-used)
##counts
paste $MY_PREVIOUS_OUTPUT_COUNTS cut_final_output_counts.txt > $MY_WORKING_DIRECTORY/edited_final_output_counts2.txt
# when I paste them oddly it adds carriage returns between the columns from the two separate files so remove these:
tr -d '\r' < edited_final_output_counts2.txt > edited_final_output_counts2_fixed.txt
mv $MY_WORKING_DIRECTORY/edited_final_output_counts2_fixed.txt $MY_WORKING_DIRECTORY/edited_final_output_counts.txt

##tpm
paste $MY_PREVIOUS_OUTPUT_TPM cut_final_output_tpm.txt > $MY_WORKING_DIRECTORY/edited_final_output_tpm2.txt
# when I paste them oddly it adds carriage returns between the columns from the two separate files so remove these:
tr -d '\r' < edited_final_output_tpm2.txt > edited_final_output_tpm2_fixed.txt
mv $MY_WORKING_DIRECTORY/edited_final_output_tpm2_fixed.txt $MY_WORKING_DIRECTORY/edited_final_output_tpm.txt




