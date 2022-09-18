rm -r ./RSA_csvs; mkdir ./RSA_csvs

for filename in ./new_templates/pMHCI_AC/*.pdb; do
	peptide_file=${filename#*/new_templates/pMHCI_AC/}
	peptide=$(cat ./template_sequences/peptide/${peptide_file} | sed -n '2 p')
	peptide_length=$(echo -n "$peptide" | wc -c)
	pdb_delchain -C ${filename} | pdb_sort | pdb_tidy | pdb_reatom | grep -E 'ATOM  |TER  |END  ' > receptor.pdb
	echo -ne "${peptide_file}\r"
	touch "./RSA_csvs/${peptide_file}"
	for ((i = 1; i <= ${peptide_length}; ++i)); do
		#echo $i
		pdb_delchain -A ${filename} | pdb_sort | pdb_tidy | pdb_reatom | grep -E 'ATOM  |TER  |END  ' | pdb_selres -${i} > residue.pdb
		pdb_merge receptor.pdb residue.pdb | pdb_sort | pdb_tidy | pdb_reatom > complex.pdb
		./naccess complex.pdb >> naccess.log
		cat complex.rsa | grep "C  " | grep "RES" | awk '{print $2, $4, $5, $6, $7, $8}' | tr ' ' ',' | awk -F, -v pep="$peptide_file", -v OFS="," '{$7=pep}1' | awk -F',' '{print $7, $1, $2, $3, $4, $5, $6}' | tr ' ' ',' >> "./RSA_csvs/${peptide_file}"
	done
done

cat "./RSA_csvs/*.pdb" > "./RSA_csvs/RSA.csv"