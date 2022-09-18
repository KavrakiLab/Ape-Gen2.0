# Delete directories when they previousy exist and remake them?
rm -r ./new_templates/allele_only; mkdir ./new_templates/allele_only
rm -r ./new_templates/peptide_only; mkdir ./new_templates/peptide_only
rm -r ./new_templates/pMHCI_A; mkdir ./new_templates/pMHCI_A
rm -r ./new_templates/pMHCI_AC; mkdir ./new_templates/pMHCI_AC

for filename in ./new_templates/pMHCI/*.pdb; do
	file=${filename#*/new_templates/pMHCI/}
	echo -ne "${file}\r"
	pdb_rplchain -M:A ${filename} > "./new_templates/pMHCI_A/${file}"
	pdb_rplchain -P:C "./new_templates/pMHCI_A/${file}" > "./new_templates/pMHCI_AC/${file}"
	length=$(cat "./new_templates/pMHCI_AC/${file}" | grep "ATOM   " | grep -E '^.{16}[A]' | wc -l)
	if [[ $length -gt 0 ]]
	then
		cat "./new_templates/pMHCI_AC/${file}" | grep -v -E '^.{16}[B]' | grep -v -E '^.{16}[C]' | grep -v -E '^.{16}[D]' | sed 's/^\(.\{16\}\)A/\1 /' | pdb_occ > "./new_templates/pMHCI_AC/${file/.pdb/_A.pdb}"
		mv "./new_templates/pMHCI_AC/${file/.pdb/_A.pdb}" "./new_templates/pMHCI_AC/${file}"	
	fi
	pdb_selchain -C "./new_templates/pMHCI_AC/${file}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "./new_templates/peptide_only/${file}"
	pdb_selchain -A "./new_templates/pMHCI_AC/${file}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "./new_templates/allele_only/${file}"
done