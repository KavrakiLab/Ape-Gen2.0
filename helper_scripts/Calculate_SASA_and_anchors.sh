#function pause(){
#   read -p "$*"
#}

#mkdir ./new_templates/SASA_csvs
#for filename in ./new_templates/pMHCI_plus_hydrogens/*.pdb; do
#	peptide_file=${filename#*/new_templates/pMHCI_plus_hydrogens/}
#	freesasa --format=rsa --radii=naccess --probe-radius 1.4 --shrake-rupley -n 200 ${filename} | grep "C  " | grep "RES" | awk '{print $4, $7, $8}' > "./new_templates/temp_pMHC.csv"
#	freesasa --format=rsa --radii=naccess --probe-radius 1.4 --shrake-rupley -n 200 "./new_templates/peptide_only/${peptide_file}" | grep "C  " | grep "RES" | awk '{print $4, $7, $8}' > "./new_templates/temp_peptide.csv"
#	join -j 1 ./new_templates/temp_pMHC.csv ./new_templates/temp_peptide.csv | tr ' ' ',' | awk -F, -v OFS="," '{$6=$4-$2}1' | awk -F, -v OFS="," '{$7=$5-$3}1' | awk -F, -v pep="$peptide_file", -v OFS="," '{$8=pep}1' | awk -F',' '{print $8, $1, $6, $7}'  > "./new_templates/temp_sidechains.csv"
#	freesasa --select="carbon_alpha, name CA" --depth=atom --format=xml --radii=naccess --probe-radius 1.4 --shrake-rupley -n 200 "./new_templates/peptide_only/${peptide_file}" | sed -n -e '/chain label="C"/,$p' | grep -Pio 'name="CA".*area="\K[^"]*' > "./new_templates/temp_peptide.csv"
#	freesasa --select="carbon_alpha, name CA" --depth=atom --format=xml --radii=naccess --probe-radius 1.4 --shrake-rupley -n 200 ${filename} | sed -n -e '/chain label="C"/,$p' | grep -Pio 'name="CA".*area="\K[^"]*' > "./new_templates/temp_pMHC.csv"
#	paste -d, ./new_templates/temp_pMHC.csv ./new_templates/temp_peptide.csv | awk -F, -v OFS="," '{$3=$2-$1}1' | awk -F',' '{print $3}' > "./new_templates/temp_CA.csv"
	#pause 'Press [Enter] key to continue...'
#	paste -d' ' "./new_templates/temp_sidechains.csv" "./new_templates/temp_CA.csv" | tr ' ' ',' | awk -F, -v OFS="," '{$3=$3+$5}1' | awk -F, -v OFS="," '{$4=$4+$5}1' | sort -t, -r -rn -k4,4 | awk -F',' '{print $1, $2, $3, $4}' > "./new_templates/SASA_csvs/${peptide_file}"
#	rm "./new_templates/temp_pMHC.csv"
#	rm "./new_templates/temp_peptide.csv"
#	rm "./new_templates/temp_CA.csv"
#done
#cat ./new_templates/SASA_csvs/*.pdb > ./new_templates/SASA_values.csv
#sed -i 1i"PDB,residue_number,Delta_SASA,Delta_rel_SASA" ./new_templates/SASA_values.csv
#rm -r ./new_templates/SASA_csvs

naccess ./new_templates/pMHCI_plus_hydrogens/1DUZ.pdb
#cat .rsa | grep "C  " | grep "RES" | awk '{print $4, $7, $8}' > "./new_templates/temp_pMHC.csv"
#rm ../Naccess/naccess2.1.1/.rsa
#naccess ./new_templates/peptide_only/1DUZ.pdb 
#cat .rsa | grep "C  " | grep "RES" | awk '{print $4, $7, $8}' > "../new_templates/temp_peptide.csv"
#rm ../Naccess/naccess2.1.1/.rsa
#join -j 1 ../new_templates/temp_pMHC.csv ../new_templates/temp_peptide.csv | tr ' ' ',' | awk -F, -v OFS="," '{$6=$4-$2}1' | awk -F, -v OFS="," '{$7=$5-$3}1' | awk -F, -v pep="$peptide_file", -v OFS="," '{$8=pep}1' | awk -F',' '{print $8, $1, $6, $7}'  > "../new_templates/temp_sidechains.csv"