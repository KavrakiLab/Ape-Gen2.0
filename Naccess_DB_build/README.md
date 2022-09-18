- prepape_templates_for_SASA_analysis_As_only.sh
Takes the PANDORA dataset and:
	1) Renames the chains
	2) Checks how many of these .pdbs have ambiguous positions (with A,B,C,Ds)
	3) and keep only the As for simplicity

- Extract_allele_and_peptide_sequences.sh:
Goes through the peptide and allele files and extracts the sequences, and puts them in .fasta files for matching (I guess)

- Calculate_RSA.sh:
For each peptide in the DB:
	For each residue in the peptide:
		Remove all other peptide residues (to mitigate the effects) and calculate RSA with the receptor present. 

- pairwise_alignment.py:
This will basically identify the allele allotypes for each allele sequence, given the MHC_seq.fasta, which is a collection of MHC sequences.
	Results are store in template_names.log

- Make_anchor_spreadsheet.py
Brings it all together and makes the DB spreadsheet, which is Template_DB_information.csv