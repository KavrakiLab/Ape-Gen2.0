import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier

from biopandas.pdb import PandasPdb
from Bio import Align

import sys
import os
import re
import pickle as pkl

from helper_scripts.Ape_gen_macros import apply_function_to_file, remove_file,  \
											rev_anchor_dictionary, all_three_to_one_letter_codes,  \
											move_file, copy_file, merge_and_tidy_pdb,              \
											replace_chains, remove_remarks_and_others_from_pdb,    \
											delete_elements, extract_CONECT_from_pdb, csp_solver,  \
											standard_three_to_one_letter_code, anchor_dictionary,  \
											verbose, extract_anchors
from classes.pMHC_class import pMHC

from subprocess import call

from pdbtools import pdb_tofasta, pdb_delelem

from openmm.app import PDBFile, ForceField, Modeller, CutoffNonPeriodic

class Peptide(object):

	def __init__(self, sequence, PTM_list=[], pdb_filename="", primary_anchors=None, secondary_anchors=None, index=-1):
		self.sequence = sequence # make sequence only the AAs
		self.PTM_list = PTM_list # have PTM list keep track of PTMs
		self.pdb_filename = pdb_filename
		self.pdbqt_filename = None
		self.primary_anchors = primary_anchors
		self.secondary_anchors = secondary_anchors
		self.index = index


	@classmethod
	def frompdb(cls, pdb_filename, primary_anchors=None, secondary_anchors=None, peptide_index=-1):
		# Initialize peptide from a .pdb file
		ppdb = PandasPdb()
		ppdb.read_pdb(pdb_filename)
		pdb_df = ppdb.df['ATOM']
		pdb_df = pdb_df[pdb_df['chain_id'] == 'C'][['residue_name', 'residue_number']].drop_duplicates()
		peptide_3letter_list = pdb_df['residue_name'].tolist()

		if len(peptide_3letter_list) == 0:
			if verbose(): print("Chain C does not exist in given .pdb file, check your format")
			
		try:
			peptide_sequence = ''.join([all_three_to_one_letter_codes[aa] for aa in peptide_3letter_list])
		except KeyError as e:
			if verbose(): print("There is something wrong with your .pdb 3-letter amino acid notation")
		return cls(pdb_filename=pdb_filename, sequence=peptide_sequence, primary_anchors=primary_anchors, secondary_anchors=secondary_anchors, index=peptide_index)

	@classmethod
	def init_peptide(cls, peptide_input):
		if verbose(): print("Processing Peptide Input")
		if peptide_input.endswith(".pdb"):
			# Fetch peptide sequence from .pdb and use that .pdb as a template --> Only when REDOCKING!
			# Maybe here have a routine that calculates the RSA? Using Naccess (Is that legal even?)
			# TODO: should anchors be set from arguments?
			#	fromPDB is only for redocking?
			#	is there never an instance where the input will only be chain C? 
			return Peptide.frompdb(peptide_input, anchors="")
		else:
			# Fetch template from peptide template list
			# peptide = Peptide.fromsequence(peptide_input)
			# peptide, template_anchors = Peptide.fromsequence(peptide_input, receptor.allotype, anchors)
			peptide_sequence_noPTM, peptide_PTM_list = PTM_processing(peptide_input)
			return cls(sequence=peptide_sequence_noPTM, PTM_list=peptide_PTM_list)
	

	def get_peptide_template(self, receptor_allotype, anchors, anchor_selection, cv=''):

		# Current policy of selecting/chosing peptide templates is:
		sequence_length = len(self.sequence)
		templates = pd.read_csv("./helper_files/Updated_template_information.csv") # Fetch template info

		if cv != '': templates = templates[~templates['pdb_code'].str.contains(cv, case=False)]
		# removes pdb code of peptide in order to cross validate (just for testing)

		# 1) Use RF to predict which anchors are to be selected (or given as an input?), and fetch the best matches
		# Let's assume for now that we have the RF, and we will be fetching templates from the DB
		# (maybe play with user input for now?)
		if anchors == "":
			
			if verbose(): print("Determining anchors for given peptide sequence and allele allotype")
			# Load the MHCflurry frequencies
			frequencies = pd.read_csv("./helper_files/mhcflurry.ba.frequency_matrices.csv")

			frequencies = frequencies[(frequencies['cutoff_fraction'] == 0.01)]
			frequencies['X'] = np.zeros(frequencies.shape[0])
			frequencies_alleles = pd.unique(frequencies['allele'])

			if receptor_allotype in frequencies_alleles:
				if verbose(): print("Receptor allotype has a known MHC binding motif!")
				anchors = extract_anchors(self.sequence, receptor_allotype, frequencies)
			else:
				print("Receptor allotype has no known MHC binding motif... Anchors are defined as canonical!")
				anchors = "2," + str(len(peptide_sequence))

		if verbose(): print("Predicted anchors for the peptide: ", anchors)
		anchors_not = process_anchors(anchors, self.sequence)
		templates['Major_anchor_not'] = templates['Major_anchor_not'].apply(lambda x: x.split(",")).apply(set) # Convert the column into a set, and do set distances
		templates['Secondary_anchor_not'] = templates['Secondary_anchor_not'].apply(lambda x: x.split(",")).apply(set) # Convert the column into a set, and do set distances
		
		# Choosing the major anchors as template choosing mechanism (secondary anchors are way to complicated)
		templates['jaccard_distance'] = templates['Major_anchor_not'].apply(lambda x: jaccard_distance(x, anchors_not))
		templates = templates[templates['jaccard_distance'] == templates['jaccard_distance'].max()].dropna()

		# 2) Bring the peptide template of MHC closer to the query one given the peptide binding motifs
		sub_alleles = pd.unique(templates['MHC']).tolist()
		sub_alleles.append("Allele")
		similarity_matrix = pd.read_csv("./helper_files/" + str(sequence_length) + "mer_similarity.csv")[sub_alleles]
		allele_of_interest = similarity_matrix[similarity_matrix["Allele"] == receptor_allotype].drop("Allele", axis=1).T
		similar_alleles = allele_of_interest[allele_of_interest == allele_of_interest.min().values[0]].dropna().index.values[0]

		templates = templates[templates['MHC'] == similar_alleles]

		# 3) Select the one that is closer in terms of anchor residues
		peptide_anchor_sequence = self.sequence[:2] + self.sequence[(sequence_length - 2):]
		template_anchor_sequences = (templates['peptide'].str[:2] + templates['peptide'].str[(sequence_length - 2):]).tolist()
		aligner = Align.PairwiseAligner()
		aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
		score_list = []
		for template_anchor_sequence in template_anchor_sequences:
			score_list.append(aligner.score(peptide_anchor_sequence, template_anchor_sequence))
		templates['anchor_score'] = score_list
		templates = templates[templates['anchor_score'] == templates['anchor_score'].max()].dropna()

		# 4) If there are duplicate results, select the one closer to the whole sequence
		score_list = []
		template_sequences = templates['peptide'].tolist()
		aligner.open_gap_score = -0.5
		aligner.extend_gap_score = -0.5
		for template_sequence in template_sequences:
			score_list.append(aligner.score(self.sequence, template_sequence))
		templates['peptide_score'] = score_list
		templates = templates[templates['peptide_score'] == templates['peptide_score'].max()].dropna()

		# 5) If there are duplicate results, select at random!
		final_selection = templates.sample(n=1)
		peptide_template_file = './new_templates/' + final_selection['pdb_code'].values[0]
		template_peptide_length = final_selection['peptide_length'].values[0]

		# 6) Extract the anchor position numbers for the anchor tolerance step!
		# CAUTION: This could cause inconsistences if the peptide sizes differ greatly, but not really, just making the anchor tolerance step a little bit more obscure
		template_major_anchors = sorted([rev_anchor_dictionary[anchor][str(template_peptide_length)] for anchor in list(final_selection['Major_anchor_not'].values[0])])
		template_secondary_anchors = sorted([rev_anchor_dictionary[anchor][str(template_peptide_length)] for anchor in list(final_selection['Secondary_anchor_not'].values[0])])
		peptide_primary_anchors = sorted([rev_anchor_dictionary[anchor][str(sequence_length)] for anchor in anchors_not])
		peptide_second_anchors = sorted([rev_anchor_dictionary[anchor][str(sequence_length)] for anchor in list(final_selection['Secondary_anchor_not'].values[0])])

		# 7) Define the peptide template object
		peptide_template = pMHC(pdb_filename=peptide_template_file, 
								peptide=Peptide.frompdb(pdb_filename=peptide_template_file, 
														  primary_anchors=template_major_anchors,
														  secondary_anchors=template_secondary_anchors))

		# 8) Return both the peptide object, as well as the peptide template that was chosen
		self.pdb_filename = None
		self.primary_anchors = peptide_primary_anchors 
		self.secondary_anchors = peptide_second_anchors
		return peptide_template

	# def add_sidechains(self, filestore):
	# 	fixer = PDBFixer(filename=self.pdb_filename)
	# 	fixer.findMissingResidues()
	# 	fixer.removeHeterogens(True) #  True keeps water molecules while removing all other heterogens, REVISIT!
	# 	fixer.findMissingAtoms()
	# 	fixer.addMissingAtoms()
	# 	fixer.addMissingHydrogens(7.0) # Ask Mauricio about those
	# 	#fixer.addSolvent(fixer.topology.getUnitCellDimensions()) # Ask Mauricio about those
	# 	self.pdb_filename = filestore + '/add_sidechains/PTMed_' + str(self.index) + '.pdb'
	# 	PDBFile.writeFile(fixer.topology, fixer.positions, open(self.pdb_filename, 'w'))

	# 	# So my hypothesis now here is that Modeller that is being used to add hydrogens, has a specification
	# 	# in its file, hydrogens.xml, that adds a methyl group of H2 and H3 (as well as the OXT) that really mess up prepare_ligard4.py.
	# 	# Let's delete those entries and see how this goes:
	# 	# UPDATE: It's probably not this, uncomment if necessary
	# 	#delete_modeller_hydrogens = delete_elements(self.pdb_filename, ["H2", "H3", "OXT"])
	# 	#overwritten = ''.join(delete_modeller_hydrogens)
	# 	#with open(self.pdb_filename, 'w') as PTMed_file:
	# 	#	PTMed_file.write(overwritten)

	# 	# Before finishing, also copy the file to the PTM floder, as the process is going to be self-referential (same input output for the PTM)
	# 	copy_file(filestore + '/add_sidechains/PTMed_' + str(self.index) + '.pdb', 
	# 						filestore + '/PTMed_peptides/PTMed_' + str(self.index) + '.pdb')

	def perform_PTM(self, filestore):
		# Unfortunately, I have to revert to stupid system calls here, because I cannot call pytms from python
		# Maybe one day...
		log_file = filestore + '/PTMed_peptides/PTM.log'
		self.pdb_filename = filestore + "/PTMed_peptides/PTMed_" + str(self.index) + ".pdb"
		for ptm in self.PTM_list:
			PTM, selection = ptm.split(' ', 1)
			call(["pymol -qc ./pymol_scripts/" + PTM + ".pml -- " + self.pdb_filename + " " + selection + " " + self.pdb_filename + " > " + log_file + " 2>&1"], shell=True)

		# For some reason, after this step, I get peptide .pdb files with:

		# A. Chain A. I want to make it into chains C as before
		apply_function_to_file(replace_chains, self.pdb_filename, chain_from="A", chain_to="C")

		# B. Weird H0 pymol hydrogens that I want to delete. This is added to other residues during the PTM, so I need to remove them
		apply_function_to_file(delete_elements, self.pdb_filename, element_set=["H0"], chains=["C"])

		# C. I need to re-organize atom indexes, which are a proper mess
		PTMed_tidied = filestore + "/PTMed_peptides/PTMed_" + str(self.index) + "tidied.pdb"
		merge_and_tidy_pdb([self.pdb_filename], PTMed_tidied)
		copy_file(PTMed_tidied, self.pdb_filename)
		remove_file(PTMed_tidied)

	def prepare_for_scoring(self, filestore, current_round, addH):
		prep_peptide_loc = "/conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
		self.pdbqt_filename = filestore + "/pdbqt_peptides/peptide_" + str(self.index) + ".pdbqt"
		clean = "lps" if addH == "all" else "nphs_lps"

		call(["python2.7 " + prep_peptide_loc + " -l " + self.pdb_filename + " -o " + self.pdbqt_filename + " -A None -Z -U " + clean + " -g -s > " + filestore + "/pdbqt_peptides/prepare_ligand4.log 2>&1"], shell=True)

		# If the resulting .pdbqt is faulty, delete it. If it does not exist, it is also faulty, so skip whatever else. 
		try:
			seq = pdb_tofasta.run(open(self.pdbqt_filename, 'r'), multi=False)
		except FileNotFoundError:
			return True
		seq = ''.join(seq).split("\n")[1]
		if(len(seq) != len(self.sequence)):
			#remove_file(self.pdbqt_filename)
			with open(filestore + "/per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(current_round) + "," + str(self.index) + ",Rejected by prepare_ligand4.py,-\n")
			return True
		else:
			return False


	def dock_score_with_SMINA(self, filestore, receptor, addH):

		# SMINA docking and scoring
		add_hydrogens = "True" if addH != "none" else "False"
		self.pdb_filename =  filestore + "/scoring_results/model_" + str(self.index) + ".pdb"
		if not receptor.useSMINA and receptor.doMinimization:
			call(["smina -q --scoring vinardo --out_flex " + filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 8 --local_only --minimize --flexres " + receptor.flexible_residues + \
				  " --energy_range 100 --addH " + add_hydrogens + " --out " + self.pdb_filename + " > " + \
				  filestore + "/scoring_results/smina.log 2>&1"], shell=True)
		elif not receptor.useSMINA and not receptor.doMinimization:
			call(["smina -q --scoring vinardo --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 8 --local_only --minimize --energy_range 100 --addH " + add_hydrogens + " --out " + self.pdb_filename + " > " + \
				  filestore + "/scoring_results/smina.log 2>&1"], shell=True)
			#move_file(receptor.pdb_filename, filestore + "/receptor_smina_min.pdb")
		elif receptor.useSMINA and receptor.doMinimization:
			call(["smina -q --out_flex " + filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 8 --local_only --minimize --flexres " + receptor.flexible_residues + \
				  " --energy_range 100 --addH " + add_hydrogens + " --out " + self.pdb_filename + " > " + \
				  filestore + "/scoring_results/smina.log 2>&1"], shell=True)
		elif receptor.useSMINA and not receptor.doMinimization:
			call(["smina -q --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 8 --local_only --minimize --energy_range 100 --addH " + add_hydrogens + " --out " + self.pdb_filename + " > " + \
				  filestore + "/scoring_results/smina.log 2>&1"], shell=True)
			#move_file(receptor.pdb_filename, filestore + "/receptor_smina_min.pdb")

	def score_with_SMINA(self, filestore, receptor):
		self.pdb_filename = filestore + "/scoring_results/model_" + str(self.index) + ".pdb" 
		call(["smina -q --score_only --ligand " + self.pdbqt_filename + \
			  " --receptor " + receptor.pdbqt_filename + " --out " + self.pdb_filename + \
			  " > " + filestore + "/scoring_results/smina.log 2>&1"], shell=True)
		move_file(receptor.pdb_filename, filestore + "/minimized_receptors/receptor_" + str(self.index) + ".pdb")    

	def compute_anchor_tolerance(self, filestore, receptor, peptide_template_anchors_xyz, anchor_tol, current_round):

		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(self.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only anchors
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(self.secondary_anchors)]
		# Only carbon-alpha atoms
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'] == 'CA']

		# Only positions
		pdb_peptide_anchors_xyz = pdb_df_peptide[['x_coord', 'y_coord', 'z_coord']].to_numpy()

		# If difference is smaller than the tolerance, keep the file, else don't
		anchor_difference = np.linalg.norm(pdb_peptide_anchors_xyz - peptide_template_anchors_xyz, axis=1)
		if np.all(anchor_difference < anchor_tol):
			
			# Keep the scores of the remaining survivors
			with open(self.pdb_filename, 'r') as peptide_handler:
				next(peptide_handler) # Skipping first line
				affinity = peptide_handler.readline().replace("\n", "").split(" ")[2]
			with open(filestore + "/per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(current_round) + "," + str(self.index) + ",Successfully Modeled," + str(affinity) + "\n")

			dst = filestore + "/anchor_filtering/peptide_" + str(self.index) + ".pdb"
			copy_file(self.pdb_filename, dst)
			self.pdb_filename = dst
			return False

		else:
			# delete the minimized receptor coming from SMINA
			if(receptor.doMinimization and os.path.exists(filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb")): remove_file(filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb")
			
			# Keep a log for the anchor difference
			with open(filestore + "/anchor_filtering/peptide_" + str(self.index) + ".log", 'a+') as anchor_log:
				anchor_log.write(str(current_round) + "," + str(self.index) + "," + ','.join(map(str, anchor_difference))) 
			
			# Keep this result for final printing
			faulty_positions = (anchor_difference > anchor_tol)*self.secondary_anchors
			faulty_positions = " and ".join(np.char.mod('%d', faulty_positions[faulty_positions != 0]))
			with open(filestore + "/per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(current_round) + "," + str(self.index) + ",Anchor tolerance violated in positions " + faulty_positions + ",-\n")
			
			return True

	def fix_flexible_residues(self, filestore, receptor, current_round, addH):

		# Make the flexible receptor output from the SMINA --out_flex argument
		#minimized_receptor_loc = filestore + "/4_SMINA_data/minimized_receptors/receptor_" + str(self.index) + ".pdb"
		#if receptor.doMinimization:
		#	call(["python ./helper_scripts/makeflex.py " + \
		#		  filestore + "/4_SMINA_data/receptor_for_smina.pdb " + \
		#		  filestore + "/4_SMINA_data/flexible_receptors/receptor_" + str(self.index) + ".pdb " + \
		#		  minimized_receptor_loc],
		#		  stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'), shell=True)

		# Alternative scenario as makeflex.py is probably unstable: Solve the CSP using the CONECT fields to determine the true identity of the atoms

		# Making the CONECT list first:
		edge_list = extract_CONECT_from_pdb(filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb")

		original_ppdb = PandasPdb()
		original_ppdb.read_pdb(filestore + "/receptor_for_smina.pdb")
		original_pdb_df = original_ppdb.df['ATOM']

		flexible_ppdb = PandasPdb()
		flexible_ppdb.read_pdb(filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb")
		flexible_pdb_df = flexible_ppdb.df['ATOM']

		# Main Routine: For each flexible residue, solve the CSP and rename the atoms based on the CONECT fields
		flexible_residues = pd.unique(flexible_pdb_df['residue_number'])
		list_of_dataframes = []
		for flex_residue in flexible_residues:
			sub_origin_pdb = original_pdb_df[(original_pdb_df['residue_number'] == flex_residue) & (original_pdb_df['chain_id'] == 'A')].copy()
			sub_pdb = flexible_pdb_df[flexible_pdb_df['residue_number'] == flex_residue].copy()
			atom_indexes = sub_pdb['atom_number'].tolist()
			sub_edge_list = [elem for elem in edge_list if ((elem[0] in atom_indexes) and (elem[1] in atom_indexes))]
			residue = (sub_pdb['residue_name'].tolist())[0]

			# Remember that CA and C are in the backbone, and are not moving at all, so their co-ordinates will be the same
			# CA location
			CA_coords = np.array(sub_origin_pdb[sub_origin_pdb['atom_name'] == 'CA'][['x_coord', 'y_coord', 'z_coord']])
			loc_indexes = np.all(np.isclose(CA_coords, np.array(sub_pdb[['x_coord', 'y_coord', 'z_coord']]), 
										  rtol=1e-05, atol=1e-08, equal_nan=False), axis=1)
			CA_loc = (sub_pdb.loc[loc_indexes, 'atom_number'].values)[0]

			# C location
			C_coords = np.array(sub_origin_pdb[sub_origin_pdb['atom_name'] == 'C'][['x_coord', 'y_coord', 'z_coord']])
			loc_indexes = np.all(np.isclose(C_coords, np.array(sub_pdb[['x_coord', 'y_coord', 'z_coord']]), 
										  rtol=1e-05, atol=1e-08, equal_nan=False), axis=1)
			C_loc = (sub_pdb.loc[loc_indexes, 'atom_number'].values)[0]

			#print(CA_loc, C_loc)
			matching = csp_solver(sub_edge_list, residue, atom_indexes, CA_loc, C_loc, addH)
			#print(matching)
			#input()
			if matching.shape[0] == 0: # Empty Solution
				# A solution was not found: Most probable case is that the CONECT fields are also broken, meaning that the conformation is invalid as it is. 
				# os.remove(filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb") ## DON'T DELETE FOR KNOW, IN CASE WE HAVE THIS ISSUE AGAIN, INSPECT THE OUTPUT
				with open(filestore + "/per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as flexible_log:
					flexible_log.write(str(current_round) + "," + str(self.index) + ",Flexible receptor conformation received was faulty,-\n") 
				return True

			sub_pdb = sub_pdb.drop(columns='atom_name').merge(matching, how='inner', on='atom_number')
			list_of_dataframes.append(sub_pdb)  

		# When done, bring the .pdb file columns to the appropriate order
		renamed_atoms = pd.concat(list_of_dataframes)
		cols = renamed_atoms.columns.tolist()
		renamed_atoms = renamed_atoms[cols[:3] + [cols[-1]] + cols[3:-1]]

		# Unify the original file with the flexible one
		flexible_ppdb.df['ATOM'] = renamed_atoms.copy()
		original_ppdb.df['ATOM'] = original_pdb_df[(~(original_pdb_df['residue_number'].isin(flexible_residues))) | (original_pdb_df['atom_name'].isin(["N", "O", "H"]))]
		original_ppdb.to_pdb(path=filestore + "/temp_" + str(self.index) + ".pdb", records=['ATOM'], gz=False, append_newline=True)
		flexible_ppdb.to_pdb(path=filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb", records=['ATOM'], gz=False, append_newline=True)
		minimized_receptor_loc = filestore + "/minimized_receptors/receptor_" + str(self.index) + ".pdb"
		merge_and_tidy_pdb([filestore + "/temp_" + str(self.index) + ".pdb", 
							filestore + "/flexible_receptors/receptor_" + str(self.index) + ".pdb"],
							minimized_receptor_loc)
		remove_file(filestore + "/temp_" + str(self.index) + ".pdb")

		return False

	def create_peptide_receptor_complexes(self, filestore, receptor):

		# Unify peptide and receptor together
		pMHC_complex = filestore + "/pMHC_complexes/pMHC_" + str(self.index) + ".pdb"
		removed = remove_remarks_and_others_from_pdb(self.pdb_filename)
		overwritten = ''.join(removed)
		with open(self.pdb_filename, 'w') as peptide_handler:
			peptide_handler.write(overwritten)
		
		if not receptor.doMinimization: copy_file(filestore + "/receptor_for_smina.pdb",
												  filestore + "/minimized_receptors/receptor_" + str(self.index) + ".pdb")
		
		merge_and_tidy_pdb([filestore + "/minimized_receptors/receptor_" + str(self.index) + ".pdb", 
							self.pdb_filename], pMHC_complex)


def AA_error_checking(amino_acid):
	if (amino_acid not in standard_three_to_one_letter_code.values()) and (amino_acid not in non_standard_three_to_one_letter_code.values()):
		if verbose(): print("The provided amino acid in the sequence is wrong")
		sys.exit(0)

def process_anchors(anchors, pep_seq):
	# Returns the set of anchors
	pep_length = len(re.sub('[a-z]', '', pep_seq)) # Remove any PTMs that may still exist in the sequence
	anchor_not = set([anchor_dictionary[str(pep_length)][str(aa_index)] for aa_index in anchors.split(",")])
	return anchor_not

def jaccard_distance(a, b):
	# Computes jaccard distance between 2 sets
	c = a.intersection(b)
	return float(len(c)) / (len(a) + len(b) - len(c))

## PTMs

# Different PTMs

phosphorylation_list = ['pS', 'pT', 'pY']
acetylation_list = ['aK'] # Check details on pytms -> This is nmot working well, it renames all hydrogens to PyMOL ones
carbamylation_list = ['cK'] # Check details on pytms
citrullination_list = ['cR'] 
methylation_list = ['mmK', 'mdK', 'mtK'] # Check details on pytms
nitration_list = ['nY', 'nW'] # Check details on pytms
s_nitrosylation_list = ['nC']
p_hydroxylation_list = ['nP'] # Check details on pytms
# malondialdehyde_list = ['maK'] # Check details on pytms, but probably it's too complicated to use that one
c_oxidation_list = ['xhC', 'xoC', 'xdC'] # Check details on pytms
m_oxidation_list = ['oM'] # Check details on pytms


def PTM_processing(sequence):
	sequence_list = re.sub( r"([A-Z])", r"\1 ", sequence).split() # Split pep sequence while retaining the PTM
		
	PTM_list = []
	sequence_noPTM = ""
	for i, amino_acid in enumerate(sequence_list):
		if len(amino_acid) > 1:
			prefix = PTM_error_checking(amino_acid)
			AA_error_checking(amino_acid[1])
			PTM_list.append(prefix + str(i + 1))
		else:
			AA_error_checking(amino_acid)
		sequence_noPTM += amino_acid[-1]
	return sequence_noPTM, PTM_list

def PTM_error_checking(amino_acid):
	prefix = amino_acid[0]
	if prefix == 'p':
		if amino_acid in phosphorylation_list:
			return "phosphorylate "
		else:
			if verbose(): print("The only amino acids that support phosphorylation are S, T and Y")
			sys.exit(0)
	elif prefix == 'n':
		if amino_acid in s_nitrosylation_list: # Keep in mind that we will need an elif here for the rest of n's
			return "nitrosylate "
		elif amino_acid in p_hydroxylation_list:
			return "proline-hydroxylate "
		elif amino_acid in nitration_list:
			return "nitrate "
		else:
			if verbose(): print("The PTMs that have the n as prefix are s-nitrosylation and p-hydroxylation (maybe nitration also). For these PTMs, the only supported amino acids are C, and P (maybe W and Y)")
			sys.exit(0)
	elif prefix == 'c':
		if amino_acid in citrullination_list: # Keep in mind that we will need an elif here for the rest of c's
			return "citrullinate "
		elif amino_acid in carbamylation_list:
			return "carbamylate "
		else:
			if verbose(): print("The PTMs that have the c as prefix are carbamylation and citrullination (maybe c_oxidation also). For these PTMs, the only supported amino acids are C, K and R")
			sys.exit(0)
	elif prefix == 'a':
		if amino_acid in acetylation_list:
			return "acetylate "
		else:
			if verbose(): print("The PTM that has the a as prefix is acetylation. For these PTM, the only supported amino acid is K.")
			sys.exit(0)
	elif prefix == 'm':
		if amino_acid in methylation_list:
			if amino_acid[1] == 'm':
				return "methylate "
			elif amino_acid[1] == 'd':
				return "di-methylate "
			elif amino_acid[1] == 't':
				return "tri-methylate "
			else:
				if verbose(): print("PTM chosen is methylation. You can only mono-methylate (mm), di-methylate (md) or tri-methylate (mt).")
				sys.exit(0)
		else:
			if verbose(): print("PTM chosen is methylation. You can only mono-methylate (mm), di-methylate (md) or tri-methylate (mt).")
			sys.exit(0)
	elif prefix == 'x':
		if amino_acid in c_oxidation_list:
			if amino_acid[1] == 'h':
				return "cysteine-hydroxydate "
			elif amino_acid[1] == 'o':
				return "cysteine-oxydate "
			elif amino_acid[1] == 'd':
				return "cysteine-dioxydate "
			else:
				if verbose(): print("PTM chosen is cysteine oxidation. You can only cysteine-hydroxidate (xh), cysteine-oxidate (xo) or cysteine-dioxidate (xd).")
				sys.exit(0)
		else:
			if verbose(): print("PTM chosen is methylation. You can only mono-methylate (mm), di-methylate (md) or tri-methylate (mt).")
			sys.exit(0)
	elif prefix == 'o':
		if amino_acid in m_oxidation_list:
			return "methionine-oxidization "
		else:
			if verbose(): print("PTM chosen is methionine oxidization. For these PTM, the only supported amino acid is M.")
			sys.exit(0)
	else:
		if verbose(): print("Wrong PTM prefix, check PTM notation")
		sys.exit(0)
