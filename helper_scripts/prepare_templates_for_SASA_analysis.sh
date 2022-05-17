for filename in ../new_templates/pMHCI/*.pdb; do
	#echo $filename
	file=${filename#*/new_templates/pMHCI/}
	echo -ne "${file}\r"
	pdb_rplchain -M:A ${filename} > "../new_templates/pMHCI_A/${file}"
	pdb_rplchain -P:C "../new_templates/pMHCI_A/${file}" > "../new_templates/pMHCI_AC/${file}"
	length=$(cat "../new_templates/pMHCI_AC/${file}" | grep "ATOM   " | grep -E '^.{16}[A]' | wc -l)
	if [[ $length -gt 0 ]]
	then
  		cat "../new_templates/pMHCI_AC/${file}" | grep -v -E '^.{16}[B]' | sed 's/^\(.\{16\}\)A/\1 /' | pdb_occ > "../new_templates/pMHCI_AC/${file/.pdb/_A.pdb}"
  		cat "../new_templates/pMHCI_AC/${file}" | grep -v -E '^.{16}[A]' | sed 's/^\(.\{16\}\)B/\1 /' | pdb_occ > "../new_templates/pMHCI_AC/${file/.pdb/_B.pdb}"
  		rm "../new_templates/pMHCI_AC/${file}"
  		pdb_selchain -C "../new_templates/pMHCI_AC/${file/.pdb/_A.pdb}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "../new_templates/peptide_only/${file/.pdb/_A.pdb}"
		pdb_selchain -C "../new_templates/pMHCI_AC/${file/.pdb/_B.pdb}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "../new_templates/peptide_only/${file/.pdb/_B.pdb}"
		pdb_selchain -A "../new_templates/pMHCI_AC/${file/.pdb/_A.pdb}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "../new_templates/anchor_only/${file/.pdb/_A.pdb}"
		pdb_selchain -A "../new_templates/pMHCI_AC/${file/.pdb/_B.pdb}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "../new_templates/anchor_only/${file/.pdb/_B.pdb}"
	else
		pdb_selchain -C "../new_templates/pMHCI_AC/${file}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "../new_templates/peptide_only/${file}"
		pdb_selchain -A "../new_templates/pMHCI_AC/${file}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "../new_templates/anchor_only/${file}"	
	fi
	#pdbfixer "./new_templates/pMHCI_AC/${file}" --output "./new_templates/pMHCI_plus_hydrogens/${file}"
	#pdb_selchain -C "./new_templates/pMHCI_plus_hydrogens/${file}" | pdb_keepcoord | pdb_reatom | pdb_tidy > "./new_templates/peptide_only/${file}"
done