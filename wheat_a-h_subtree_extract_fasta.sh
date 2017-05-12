cat NAC-a_wheat_barley_rice.txt | grep TRIAE | sed s/"TRIAE CS42 .* TGACv1 ...... AA"//g  > NAC-a_wheat.txt 
 grep -A1 -f NAC-a_wheat.txt ../wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-a_wheat.fa
 
 cat NAC-b_wheat_barley_rice.txt | grep TRIAE | sed s/"TRIAE CS42 .* TGACv1 ...... AA"//g  > NAC-b_wheat.txt 
 grep -A1 -f NAC-b_wheat.txt ../wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-b_wheat.fa
 
 cat NAC-c_wheat_barley_rice.txt | grep TRIAE | sed s/"TRIAE CS42 .* TGACv1 ...... AA"//g  > NAC-c_wheat.txt 
 grep -A1 -f NAC-c_wheat.txt ../wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-c_wheat.fa
 
 cat NAC-d_wheat_barley_rice.txt | grep TRIAE | sed s/"TRIAE CS42 .* TGACv1 ...... AA"//g  > NAC-d_wheat.txt 
 grep -A1 -f NAC-d_wheat.txt ../wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-d_wheat.fa
 
 cat NAC-e_wheat_barley_rice.txt | grep TRIAE | sed s/"TRIAE CS42 .* TGACv1 ...... AA"//g  > NAC-e_wheat.txt 
 grep -A1 -f NAC-e_wheat.txt ../wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-e_wheat.fa
 
 cat NAC-f_wheat_barley_rice.txt | grep TRIAE | sed s/"TRIAE CS42 .* TGACv1 ...... AA"//g  > NAC-f_wheat.txt 
 grep -A1 -f NAC-f_wheat.txt ../wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-f_wheat.fa
 
 cat NAC-g_wheat_barley_rice.txt | grep TRIAE | sed s/"TRIAE CS42 .* TGACv1 ...... AA"//g  > NAC-g_wheat.txt 
 grep -A1 -f NAC-g_wheat.txt ../wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-g_wheat.fa
 
 cat NAC-h_wheat_barley_rice.txt | grep TRIAE | sed s/"TRIAE CS42 .* TGACv1 ...... AA"//g  > NAC-h_wheat.txt 
 grep -A1 -f NAC-h_wheat.txt ../wheat_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-h_wheat.fa
 

 
 # now do the same for wheat, barley, rice
 
cat NAC-a_wheat_barley_rice.txt | sed s/"TRIAE CS42 .* TGACv1 ...... AA//g" | sed s/"LOC Os[0-1]"//g | sed s/MLO//g | sed s/" "/_/g > NAC-a_wheat_barley_rice_shortID.txt
grep -A1 -f NAC-a_wheat_barley_rice_shortID.txt ../wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-a_wheat_barley_rice.fa

cat NAC-b_wheat_barley_rice.txt | sed s/"TRIAE CS42 .* TGACv1 ...... AA//g" | sed s/"LOC Os[0-1]"//g | sed s/MLO//g | sed s/" "/_/g > NAC-b_wheat_barley_rice_shortID.txt
grep -A1 -f NAC-b_wheat_barley_rice_shortID.txt ../wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-b_wheat_barley_rice.fa

cat NAC-c_wheat_barley_rice.txt | sed s/"TRIAE CS42 .* TGACv1 ...... AA//g" | sed s/"LOC Os[0-1]"//g | sed s/MLO//g | sed s/" "/_/g > NAC-c_wheat_barley_rice_shortID.txt
grep -A1 -f NAC-c_wheat_barley_rice_shortID.txt ../wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-c_wheat_barley_rice.fa

cat NAC-d_wheat_barley_rice.txt | sed s/"TRIAE CS42 .* TGACv1 ...... AA//g" | sed s/"LOC Os[0-1]"//g | sed s/MLO//g | sed s/" "/_/g > NAC-d_wheat_barley_rice_shortID.txt
grep -A1 -f NAC-d_wheat_barley_rice_shortID.txt ../wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-d_wheat_barley_rice.fa

cat NAC-e_wheat_barley_rice.txt | sed s/"TRIAE CS42 .* TGACv1 ...... AA//g" | sed s/"LOC Os[0-1]"//g | sed s/MLO//g | sed s/" "/_/g > NAC-e_wheat_barley_rice_shortID.txt
grep -A1 -f NAC-e_wheat_barley_rice_shortID.txt ../wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-e_wheat_barley_rice.fa

cat NAC-f_wheat_barley_rice.txt | sed s/"TRIAE CS42 .* TGACv1 ...... AA//g" | sed s/"LOC Os[0-1]"//g | sed s/MLO//g | sed s/" "/_/g > NAC-f_wheat_barley_rice_shortID.txt
grep -A1 -f NAC-f_wheat_barley_rice_shortID.txt ../wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-f_wheat_barley_rice.fa

cat NAC-g_wheat_barley_rice.txt | sed s/"TRIAE CS42 .* TGACv1 ...... AA//g" | sed s/"LOC Os[0-1]"//g | sed s/MLO//g | sed s/" "/_/g > NAC-g_wheat_barley_rice_shortID.txt
grep -A1 -f NAC-g_wheat_barley_rice_shortID.txt ../wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-g_wheat_barley_rice.fa

cat NAC-h_wheat_barley_rice.txt | sed s/"TRIAE CS42 .* TGACv1 ...... AA//g" | sed s/"LOC Os[0-1]"//g | sed s/MLO//g | sed s/" "/_/g > NAC-h_wheat_barley_rice_shortID.txt
grep -A1 -f NAC-h_wheat_barley_rice_shortID.txt ../wheat_barley_rice_NAC_msa_NAC_domain_no_numbers_10perc_shortID.fa > NAC-h_wheat_barley_rice.fa

