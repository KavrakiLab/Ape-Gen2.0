function pause(){
   read -p "$*"
}

rm -r ./template_sequences/
mkdir ./template_sequences/
mkdir ./template_sequences/allele
mkdir ./template_sequences/peptide

for filename in ./new_templates/allele_only/*.pdb; do
   allele_file=${filename#*/new_templates/allele_only/}
   echo ">${allele_file}" > "./template_sequences/allele/${allele_file}"
   bash extract_sequence_from_pdb.sh ${filename} >> "./template_sequences/allele/${allele_file}"
done

for filename in ./new_templates/peptide_only/*.pdb; do
   peptide_file=${filename#*/new_templates/peptide_only/}
   echo ">${peptide_file}" > "./template_sequences/peptide/${peptide_file}"
   bash extract_sequence_from_pdb.sh ${filename} >> "./template_sequences/peptide/${peptide_file}"
done