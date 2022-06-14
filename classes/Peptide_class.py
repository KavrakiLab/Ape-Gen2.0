import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier

from biopandas.pdb import PandasPdb
from Bio import Align

import sys
import os
import re
import pickle as pkl

from helper_scripts.Ape_gen_macros import extract_features, rev_anchor_dictionary, all_three_to_one_letter_codes, move_file, copy_file, merge_and_tidy_pdb, replace_chains, remove_remarks_and_others_from_pdb, delete_elements, extract_CONECT_from_pdb, csp_solver, process_anchors, jaccard_distance

from subprocess import call

from pdbtools import pdb_tofasta, pdb_delelem

from openmm.app import PDBFile, ForceField, Modeller, CutoffNonPeriodic

class Peptide(object):

	def __init__(self, pdb_filename, sequence, anchors):
		self.sequence = sequence
		self.pdb_filename = pdb_filename
		self.pdbqt_filename = None
		self.anchors = anchors

	@classmethod
	def frompdb(cls, pdb_filename, anchors):
		# Initialize peptide from a .pdb file
		ppdb = PandasPdb()
		ppdb.read_pdb(pdb_filename)
		pdb_df = ppdb.df['ATOM']
		pdb_df = pdb_df[pdb_df['chain_id'] == 'C'][['residue_name', 'residue_number']].drop_duplicates()
		peptide_3letter_list = pdb_df['residue_name'].tolist()

		if len(peptide_3letter_list) == 0:
			print("Chain C does not exist in given .pdb file, check your format")
			
		try:
			peptide_sequence = ''.join([all_three_to_one_letter_codes[aa] for aa in peptide_3letter_list])
		except KeyError as e:
			print("There is something wrong with your .pdb 3-letter amino acid notation")
		return cls(pdb_filename = pdb_filename, sequence = peptide_sequence, anchors = anchors)
	
	@classmethod
	def fromsequence(cls, peptide_sequence):
		# Initialize peptide from a sequence -> Fetch template!
		templates = pd.read_csv("./template_files/n-mer-templates.csv")
		sequence_length = len(re.sub('[a-z]', '', peptide_sequence)) # Remove PTMs when fetching the template
		peptide_template = templates[templates['Pep_Length'] == sequence_length]['Template'].values[0]
		return cls(pdb_filename = ('./templates/' + peptide_template), sequence = peptide_sequence)

	@classmethod
	def fromsequence(cls, peptide_sequence, receptor_allotype, anchors):

		# Current policy of selecting/chosing peptide templates is:
		peptide_sequence_noPTM = re.sub('[a-z]', '', peptide_sequence) # Remove PTMs when fetching the template
		sequence_length = len(re.sub('[a-z]', '', peptide_sequence_noPTM)) 
		templates = pd.read_csv("./helper_files/Template_Information_notation.csv") # Fetch template info

		# 1) Use RF to predict which anchors are to be selected (or given as an input?), and fetch the best matches
		# Let's assume for now that we have the RF, and we will be fetching templates from the DB
		# (maybe play with user input for now?)
		if anchors == "":
			
			print("Determining anchors for given peptide sequence and allele allotype")
			# Load the MHCflurry frequencies
			frequencies = pd.read_csv("./helper_files/mhcflurry.ba.frequency_matrices.csv")
			frequencies = frequencies[(frequencies['cutoff_fraction'] == 0.01)]
			frequencies['X'] = np.zeros(frequencies.shape[0])
			frequencies_alleles = pd.unique(frequencies['allele'])

			if receptor_allotype in frequencies_alleles:
				print("Receptor allotype has a known MHC binding motif!")
				peptide_features = extract_features(peptide_sequence_noPTM, receptor_allotype, frequencies)
				anchor_predictors = pkl.load(open("./helper_files/anchor_predictors.pkl", "rb"))
				predictions = np.zeros(sequence_length)
				for i in range(len(anchor_predictors[0])):
					predictions += anchor_predictors[0][i].predict(np.vstack(peptide_features))
				predictions /= len(anchor_predictors[0])
				anchors = list(np.argwhere(predictions >= 0.5).flatten() + 1)
				anchors = ",".join(map(str, anchors))
			else:
				print("Receptor allotype has no known MHC binding motif... Anchors are defined as canonical!")
				anchors = "2,9"

		print("Predicted anchors for the peptide: ", anchors)
		anchors_not = process_anchors(anchors, peptide_sequence_noPTM)
		templates['anchor_not'] = templates['anchor_not'].apply(lambda x: x.split(",")).apply(set) # Convert the column into a set, and do set distances
		templates['jaccard_distance'] = templates['anchor_not'].apply(lambda x: jaccard_distance(x, anchors_not))
		templates = templates[templates['jaccard_distance'] == templates['jaccard_distance'].max()].dropna()

		# 2) Bring the peptide template of MHC closer to the query one given the peptide binding motifs
		sub_alleles = pd.unique(templates['MHC']).tolist()
		sub_alleles.append("Allele")
		similarity_matrix = pd.read_csv("./helper_files/" + str(sequence_length) + "mer_similarity.csv")[sub_alleles]
		allele_of_interest = similarity_matrix[similarity_matrix["Allele"] == receptor_allotype].drop("Allele", axis=1).T
		similar_alleles = allele_of_interest[allele_of_interest == allele_of_interest.min().values[0]].dropna().index.values[0]
		templates = templates[templates['MHC'] == similar_alleles]

		# 3) Select the one that is closer in terms of anchor residues
		peptide_anchor_sequence = peptide_sequence_noPTM[:2] + peptide_sequence_noPTM[(sequence_length - 2):]
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
		for template_sequence in template_sequences:
			score_list.append(aligner.score(peptide_sequence_noPTM, template_sequence))
		templates['peptide_score'] = score_list
		templates = templates[templates['peptide_score'] == templates['peptide_score'].max()].dropna()

		# 5) If there are duplicate results, select at random!
		final_selection = templates.sample(n=1)
		peptide_template = final_selection['pdb_code'].values[0]
		template_peptide_length = final_selection['peptide_length'].values[0]
		template_anchors_not = final_selection['anchor_not'].values[0] 
		
		# 6) Before the end, it is a good idea here to match the predicted/set anchors of the peptide
		# with the anchors of the template for the anchor tolerance step. General workflow is:
		# A) Take intersection of anchor notation (this will force equal number of anchors)
		# B) Sort
		# C) Extract the numbers using the peptide lengths
		anchor_union = list(anchors_not.intersection(template_anchors_not))
		template_anchors = sorted([rev_anchor_dictionary[anchor][str(template_peptide_length)] for anchor in anchor_union])
		peptide_anchors = sorted([rev_anchor_dictionary[anchor][str(sequence_length)] for anchor in anchor_union])

		return cls(pdb_filename = ('./new_templates/' + peptide_template), 
				   sequence = peptide_sequence, 
				   anchors = peptide_anchors), template_anchors


	def perform_PTM(self, filestore, peptide_index, PTM_list):
		# Unfortunately, I have to revert to stupid system calls here, because I cannot call pytms from python
		# Maybe one day...
		log_file = filestore + '/PTMed_peptides/PTM.log'
		self.pdb_filename = filestore + "/PTMed_peptides/PTMed_" + str(peptide_index) + ".pdb"
		for ptm in PTM_list:
			PTM, selection = ptm.split(' ', 1)
			call(["pymol -qc ./pymol_scripts/" + PTM + ".pml -- " + self.pdb_filename + " " + selection + " " + self.pdb_filename + " > " + log_file + " 2>&1"], shell=True)

		# For some reason, after this step, I get peptide .pdb files with:

		# A. Chain A. I want to make it into chains C as before
		rechained = replace_chains(self.pdb_filename, "A", "C")
		overwritten = ''.join(rechained)
		with open(self.pdb_filename, 'w') as PTMed_file:
			PTMed_file.write(overwritten)

		# B. A weird H01 pymol hydrogen that I want to delete. This is added to other residues during the PTM, so I need to remove it:
		residue_list = list(range(1, len(self.sequence)))
		for ptm in PTM_list:
			PTM, selection = ptm.split(' ', 1)
			residue_list.remove(int(selection))

		delete_pymol_residues = delete_elements(self.pdb_filename, ["H01"], chains=["C"], residues=residue_list)
		overwritten_2 = ''.join(delete_pymol_residues)
		with open(self.pdb_filename, 'w') as PTMed_file:
			PTMed_file.write(overwritten_2)

		# C. I need to re-organize atom indexes, which are a proper mess
		PTMed_tidied = filestore + "/PTMed_peptides/PTMed_" + str(peptide_index) + "tidied.pdb"
		merge_and_tidy_pdb([self.pdb_filename], PTMed_tidied)
		copy_file(PTMed_tidied, self.pdb_filename)
		os.remove(PTMed_tidied)

	def prepare_for_scoring(self, filestore, peptide_index, current_round):
		prep_peptide_loc = "/conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
		self.pdbqt_filename = filestore + "/pdbqt_peptides/peptide_" + str(peptide_index) + ".pdbqt"
		call(["python2.7 " + prep_peptide_loc + " -l " + self.pdb_filename + " -o " + self.pdbqt_filename + " -A None -Z -U lps -g -s > " + filestore + "/pdbqt_peptides/prepare_ligand4.log 2>&1"], shell=True)

		# If the resulting .pdbqt is faulty, delete it. If it does not exist, it is also faulty, so skip whatever else. 
		try:
			seq = pdb_tofasta.run(open(self.pdbqt_filename, 'r'), multi=False)
		except FileNotFoundError:
			return True
		seq = ''.join(seq).split("\n")[1]
		if(len(seq) != len(self.sequence)):
			#os.remove(self.pdbqt_filename)
			with open(filestore + "/per_peptide_results/peptide_" + str(peptide_index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(current_round) + "," + str(peptide_index) + ",Rejected by prepare_ligand4.py,-\n")
			return True
		else:
			return False

	def dock_score_with_SMINA(self, filestore, receptor, peptide_index):

		# SMINA docking and scoring
		self.pdb_filename =  filestore + "/Scoring_results/model_" + str(peptide_index) + ".pdb"
		if not receptor.useSMINA and receptor.doMinimization:
			call(["smina -q --scoring vinardo --out_flex " + filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 8 --local_only --minimize --flexres " + receptor.flexible_residues + \
				  " --energy_range 100 --out " + self.pdb_filename + " > " + \
				  filestore + "/Scoring_results/smina.log 2>&1"], shell=True)
		elif not receptor.useSMINA and not receptor.doMinimization:
			call(["smina -q --scoring vinardo --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 8 --local_only --minimize --energy_range 100 --out " + self.pdb_filename + " > " + \
				  filestore + "/Scoring_results/smina.log 2>&1"], shell=True)
			#move_file(receptor.pdb_filename, filestore + "/receptor_smina_min.pdb")
		elif receptor.useSMINA and receptor.doMinimization:
			call(["smina -q --out_flex " + filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 8 --local_only --minimize --flexres " + receptor.flexible_residues + \
				  " --energy_range 100 --out " + self.pdb_filename + " > " + \
				  filestore + "/Scoring_results/smina.log 2>&1"], shell=True)
		elif receptor.useSMINA and not receptor.doMinimization:
			call(["smina -q --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 8 --local_only --minimize --energy_range 100 --out " + self.pdb_filename + " > " + \
				  filestore + "/Scoring_results/smina.log 2>&1"], shell=True)
			#move_file(receptor.pdb_filename, filestore + "/receptor_smina_min.pdb")

	def score_with_SMINA(self, filestore, receptor, peptide_index):
		self.pdb_filename = filestore + "/Scoring_results/model_" + str(peptide_index) + ".pdb" 
		call(["smina -q --score_only --ligand " + self.pdbqt_filename + \
			  " --receptor " + receptor.pdbqt_filename + " --out " + self.pdb_filename + \
			  " > " + filestore + "/Scoring_results/smina.log 2>&1"], shell=True)
		move_file(receptor.pdb_filename, filestore + "/minimized_receptors/receptor_" + str(peptide_index) + ".pdb")    

	def compute_anchor_tolerance(self, filestore, receptor, peptide_template_anchors_xyz, anchor_tol, peptide_index, current_round):

		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(self.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only anchors
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(self.anchors)]
		
		# Only carbon-alpha atoms
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'] == 'CA']
		
		# Only positions
		pdb_peptide_anchors_xyz = pdb_df_peptide[['x_coord', 'y_coord', 'z_coord']].to_numpy()

		# If difference is smaller than the tolerance, keep the file, else don't
		anchor_difference = np.linalg.norm(pdb_peptide_anchors_xyz - peptide_template_anchors_xyz, axis = 1)
		if np.all(anchor_difference < anchor_tol):
			
			# Keep the scores of the remaining survivors
			with open(self.pdb_filename, 'r') as peptide_handler:
				next(peptide_handler) # Skipping first line
				affinity = peptide_handler.readline().replace("\n", "").split(" ")[2]
			with open(filestore + "/per_peptide_results/peptide_" + str(peptide_index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(current_round) + "," + str(peptide_index) + ",Successfully Modeled," + str(affinity) + "\n")

			dst = filestore + "/Anchor_filtering/peptide_" + str(peptide_index) + ".pdb"
			copy_file(self.pdb_filename, dst)
			self.pdb_filename = dst
			return False

		else:
			# delete the minimized receptor coming from SMINA
			if(receptor.doMinimization and os.path.exists(filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb")): os.remove(filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb")
			
			# Keep a log for the anchor difference
			with open(filestore + "/Anchor_filtering/peptide_" + str(peptide_index) + ".log", 'a+') as anchor_log:
				anchor_log.write(str(current_round) + "," + str(peptide_index) + "," + ','.join(map(str, anchor_difference))) 
			
			# Keep this result for final printing
			faulty_positions = (anchor_difference > anchor_tol)*self.anchors
			faulty_positions = " and ".join(np.char.mod('%d', faulty_positions[faulty_positions != 0]))
			with open(filestore + "/per_peptide_results/peptide_" + str(peptide_index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(current_round) + "," + str(peptide_index) + ",Anchor tolerance violated in positions " + faulty_positions + ",-\n")
			
			return True

	def fix_flexible_residues(self, filestore, receptor, peptide_index, current_round):

		# Make the flexible receptor output from the SMINA --out_flex argument
		#minimized_receptor_loc = filestore + "/SMINA_data/minimized_receptors/receptor_" + str(peptide_index) + ".pdb"
		#if receptor.doMinimization:
		#   call(["python ./helper_scripts/makeflex.py " + \
		#         filestore + "/SMINA_data/receptor_for_smina.pdb " + \
		#         filestore + "/SMINA_data/flexible_receptors/receptor_" + str(peptide_index) + ".pdb " + \
		#         minimized_receptor_loc],
		#         stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'), shell=True)

		# Alternative scenario as makeflex.py is probably unstable: Solve the CSP using the CONECT fields to determine the true identity of the atoms

		# Making the CONECT list first:
		edge_list = extract_CONECT_from_pdb(filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb")

		original_ppdb = PandasPdb()
		original_ppdb.read_pdb(filestore + "/receptor_for_smina.pdb")
		original_pdb_df = original_ppdb.df['ATOM']

		flexible_ppdb = PandasPdb()
		flexible_ppdb.read_pdb(filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb")
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
										  rtol=1e-05, atol=1e-08, equal_nan=False), axis = 1)
			CA_loc = (sub_pdb.loc[loc_indexes, 'atom_number'].values)[0]

			# C location
			C_coords = np.array(sub_origin_pdb[sub_origin_pdb['atom_name'] == 'C'][['x_coord', 'y_coord', 'z_coord']])
			loc_indexes = np.all(np.isclose(C_coords, np.array(sub_pdb[['x_coord', 'y_coord', 'z_coord']]), 
										  rtol=1e-05, atol=1e-08, equal_nan=False), axis = 1)
			C_loc = (sub_pdb.loc[loc_indexes, 'atom_number'].values)[0]

			#print(CA_loc, C_loc)
			matching = csp_solver(sub_edge_list, residue, atom_indexes, CA_loc, C_loc)
			#print(matching)
			#input()
			if matching.shape[0] == 0: # Empty Solution
				# A solution was not found: Most probable case is that the CONECT fields are also broken, meaning that the conformation is invalid as it is. 
				os.remove(filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb")
				with open(filestore + "/per_peptide_results/peptide_" + str(peptide_index) + ".log", 'a+') as flexible_log:
					flexible_log.write(str(current_round) + "," + str(peptide_index) + ",Flexible receptor conformation received was faulty,-\n") 
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
		original_ppdb.to_pdb(path=filestore + "/temp_" + str(peptide_index) + ".pdb", records=['ATOM'], gz=False, append_newline=True)
		flexible_ppdb.to_pdb(path=filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb", records=['ATOM'], gz=False, append_newline=True)
		minimized_receptor_loc = filestore + "/minimized_receptors/receptor_" + str(peptide_index) + ".pdb"
		merge_and_tidy_pdb([filestore + "/temp_" + str(peptide_index) + ".pdb", 
							filestore + "/flexible_receptors/receptor_" + str(peptide_index) + ".pdb"],
							minimized_receptor_loc)
		os.remove(filestore + "/temp_" + str(peptide_index) + ".pdb")

		return False

	def create_peptide_receptor_complexes(self, filestore, receptor, peptide_index):

		# Unify peptide and receptor together
		pMHC_complex = filestore + "/pMHC_complexes/pMHC_" + str(peptide_index) + ".pdb"
		removed = remove_remarks_and_others_from_pdb(self.pdb_filename)
		overwritten = ''.join(removed)
		with open(self.pdb_filename, 'w') as peptide_handler:
			peptide_handler.write(overwritten)
		
		if not receptor.doMinimization: copy_file(filestore + "/receptor_for_smina.pdb",
												  filestore + "/minimized_receptors/receptor_" + str(peptide_index) + ".pdb")
		
		merge_and_tidy_pdb([filestore + "/minimized_receptors/receptor_" + str(peptide_index) + ".pdb", 
							self.pdb_filename], pMHC_complex)