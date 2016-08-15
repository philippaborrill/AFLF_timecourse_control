# LSBATCH: User input
#!/bin/bash
#
# LSF batch script to sort and then combine results files # 
#BSUB -q NBI-Prod128 
#BSUB -J sort_and_merge 
#BSUB -R "rusage[mem=100000]"


##############
# usage: make sure you have the previous datasets in /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples (keep a back-up copy somewhere else too!)
# these should be files: edited_final_output_counts.txt edited_final_output_tpm.txt
# make sure new samples to add (folders containing kallisto results) are in directory: /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples

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

##### move the results dir e.g. SRR1958738 to a new dir called add_samples

# move to top directory  
cd /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples

# need file with list of gene names in 1 column - in sorted order - as defined by unix with header: transcript
MY_TRANSCRIPTS_ORDERED=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/transcripts_ordered.txt

# make empty output + temp files:
echo "" > /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/output_counts.txt
echo "" > /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/output_tpm.txt
echo "" > /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/temp.txt
echo "" > /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/col4.txt
echo "" > /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/col5.txt

#  three output files: output_counts.txt, output_fpkm.txt and output_tpm.txt
MY_OUTPUT_COUNTS=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/output_counts.txt
MY_OUTPUT_TPM=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/output_tpm.txt

#temp file:
MY_TEMP_FILE=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/temp.txt

#col8, 11 and 15 files:
MY_COL4=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/col4.txt
MY_COL5=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/add_samples/col5.txt


# put the list of ordered gene names as the left hand column in the output files:
cp $MY_TRANSCRIPTS_ORDERED $MY_OUTPUT_COUNTS
cp $MY_TRANSCRIPTS_ORDERED $MY_OUTPUT_TPM

# loop through each dir ($dir_name) in turn which has a results file in it and put the values of counts and tpm into relevant output files

	for dir_name in */ ; do
	#for dir_name in `ls -d */`; do
	
   	echo "directory name is: $dir_name \n"
	
	
	cd $dir_name
	# check results.xprs has 103,275 lines
	MY_WC=`cat abundance.tsv | wc -l`
	echo "My line count is: $MY_WC \n"
	
		
		if [ "$MY_WC" == "103275" ]

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

######## don't need this currently so write exit ## can use if it want to add extra samples to files

exit

# previous samples are stored in:
MY_WORKING_DIR=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/

# have 2 files there:
MY_PREVIOUS_OUTPUT_COUNTS=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/edited_final_output_counts.txt
MY_PREVIOUS_OUTPUT_TPM=/nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/edited_final_output_tpm.txt

# remove first column from new data sheets
cut -f 2- final_output_counts.txt > cut_final_output_counts.txt
cut -f 2- final_output_tpm.txt > cut_final_output_tpm.txt

# combine the new and old data sheets in the folder "working_dir" and renamed as edited_final_output...txt (so this script can be re-used)
##counts
paste $MY_PREVIOUS_OUTPUT_COUNTS cut_final_output_counts.txt > /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/edited_final_output_counts2.txt
# when I paste them oddly it adds carriage returns between the columns from the two separate files so remove these:
tr -d '\r' < edited_final_output_counts2.txt > edited_final_output_counts2_fixed.txt
mv /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/edited_final_output_counts2_fixed.txt /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/edited_final_output_counts.txt

##tpm
paste $MY_PREVIOUS_OUTPUT_TPM cut_final_output_tpm.txt > /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/edited_final_output_tpm2.txt
# when I paste them oddly it adds carriage returns between the columns from the two separate files so remove these:
tr -d '\r' < edited_final_output_tpm2.txt > edited_final_output_tpm2_fixed.txt
mv /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/edited_final_output_tpm2_fixed.txt /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/PB_AFLF/control_timecourse/kallisto_analysis/edited_final_output_tpm.txt




