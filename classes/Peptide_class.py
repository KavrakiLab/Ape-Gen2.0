import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from math import copysign

from biopandas.pdb import PandasPdb

from Bio import Align
from Bio import pairwise2

import sys
import os
import re
import pickle as pkl
from copy import deepcopy

from helper_scripts.Ape_gen_macros import apply_function_to_file, remove_file,                     \
											rev_anchor_dictionary, all_three_to_one_letter_codes,  \
											move_file, copy_file, merge_and_tidy_pdb,              \
											replace_chains, remove_remarks_and_others_from_pdb,    \
											delete_elements, extract_CONECT_from_pdb, csp_solver,  \
											standard_three_to_one_letter_code, anchor_dictionary,  \
											verbose, count_number_of_atoms, score_sequences,       \
											predict_anchors_PMBEC, anchor_alignment,               \
											calculate_anchors_given_alignment 

from classes.pMHC_class import pMHC

from subprocess import call

from pdbtools import pdb_tofasta, pdb_delelem

from openmm.app import PDBFile, ForceField, Modeller, CutoffNonPeriodic

class Peptide(object):

	def __init__(self, sequence, PTM_list=[], pdb_filename="", primary_anchors=None, secondary_anchors=None, 
				 tilted_sequence=None, index=-1):
		self.sequence = sequence # make sequence only the AAs
		self.PTM_list = PTM_list # have PTM list keep track of PTMs
		self.pdb_filename = pdb_filename
		self.pdbqt_filename = None
		self.primary_anchors = primary_anchors
		self.secondary_anchors = secondary_anchors
		self.tilted_sequence = tilted_sequence
		self.index = index

	@classmethod
	def init_peptide(cls, peptide_input):
		if peptide_input.endswith(".pdb"):
			# Fetch peptide sequence from .pdb and use that .pdb as a template --> Only when REDOCKING!
			# Maybe here have a routine that calculates the RSA? Using Naccess (Is that legal even?)
			# TODO: should anchors be set from arguments?
			#	fromPDB is only for redocking?
			#	is there never an instance where the input will only be chain C? 
			return Peptide.frompdb(peptide_input, anchors="")
		else:
			# Fetch template from peptide template list
			peptide_sequence_noPTM, peptide_PTM_list = PTM_processing(peptide_input)
			return cls(sequence=peptide_sequence_noPTM, PTM_list=peptide_PTM_list)
	
	@classmethod
	def frompdb(cls, pdb_filename, PTM_list=[], primary_anchors=None, secondary_anchors=None, 
		        tilted_sequence=None, peptide_index=-1):
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
		return cls(pdb_filename=pdb_filename, sequence=peptide_sequence, PTM_list=PTM_list, 
			       primary_anchors=primary_anchors, secondary_anchors=secondary_anchors, 
			       tilted_sequence=tilted_sequence, index=peptide_index)

	def get_peptide_templates(self, receptor_allotype, anchors, max_no_templates, similarity_threshold, cv=''):

		if verbose(): print("\nProcessing Peptide Input: " + self.sequence)

		sequence_length = len(self.sequence)
		templates = pd.read_csv("./helper_files/Proper_files/Template_DB.csv") # Fetch template info

		# removes pdb code of peptide in order to cross validate (just for testing)
		if cv != '': templates = templates[~templates['pdb_code'].str.contains(cv, case=False)]

		# Current policy of selecting/chosing peptide templates is:
		# 1) Feature filtering to predict which are the anchors (when they are not given)
		if anchors == "":	
			if verbose(): print("Determining anchors for given peptide sequence and allele allotype")
			anchors, anchor_status = predict_anchors_PMBEC(self.sequence, receptor_allotype)
			if verbose():
				if anchor_status == "Not Known":
					print("SMM matrices could NOT be used... Defaulting to peptide sequence alignment")
				else:
					print("Receptor allotype has a SMM matrix!")

		if verbose() and anchor_status == "Known": print("Predicted anchors for the peptide: ", anchors)
		
		# Bring the templates having same anchor distance as the give peptide-MHC (NOTE: Calculate this offline as an improvement?)
		int_anchors = [int(pos) for pos in anchors.split(",")]
		diff = abs(int_anchors[0] - int_anchors[1]) 
		templates['anchor_diff'] = templates['Major_anchors'].apply(lambda x: x.split(","))
		templates['anchor_diff'] = abs(templates['anchor_diff'].str[0].astype(int) - templates['anchor_diff'].str[1].astype(int))
		if anchor_status == "Known":
			templates = templates[templates['anchor_diff'] == diff]

		# 2) Bring the peptide template of MHC closer to the query one given the peptide binding motifs + peptide similarity
		# This helps mitigating the effect of overfitting to a peptide sequence or an allele

		# 2a. Allele similarity:
		sub_alleles = pd.unique(templates['MHC']).tolist()
		sub_alleles.append("Allele")
		similarity_matrix = pd.read_csv("./helper_files/" + str(max(8, sequence_length)) + "mer_similarity.csv")[sub_alleles]
		similarity_matrix = similarity_matrix[similarity_matrix["Allele"] == receptor_allotype].drop("Allele", axis=1).T.reset_index(drop=False)
		similarity_matrix.columns = ['MHC', 'MHC_similarity']
		templates = templates.merge(similarity_matrix, how='inner', on='MHC')

		templates['Major_anchor_1'] = templates['Major_anchors'].apply(lambda x: x.split(",")).str[0].astype(int)
		templates['Major_anchor_2'] = templates['Major_anchors'].apply(lambda x: x.split(",")).str[1].astype(int)

		if anchor_status == "Known":
			templates['Anchor_diff_1'] = templates['Major_anchor_1'] - int_anchors[0]
			templates['Anchor_diff_2'] = templates['Major_anchor_2'] - templates['peptide_length'] + int_anchors[1] - sequence_length

		# 2b. Peptide similarity between anchors:
		score_list = []
		template_sequences = templates['peptide'].tolist()
		anchor_1_list = templates['Major_anchor_1'].tolist()
		anchor_2_list = templates['Major_anchor_2'].tolist()
		blosum_62 = Align.substitution_matrices.load("BLOSUM62")
		self_score = score_sequences(self.sequence, self.sequence, 0, 0, matrix=blosum_62, gap_penalty=0, norm=1) # Self score that does not care about anchor placements.
		if anchor_status == "Known":
			Anchor_diff_1 = templates['Anchor_diff_1'].tolist()
			Anchor_diff_2 = templates['Anchor_diff_2'].tolist()
			for i, template_sequence in enumerate(template_sequences):
				temp_sequence_in_question, temp_template_sequence = anchor_alignment(self.sequence, template_sequence, 
																				 Anchor_diff_1[i], Anchor_diff_2[i])
				score_list.append(score_sequences(temp_sequence_in_question, temp_template_sequence,
												  anchor_1_list[i], anchor_2_list[i],  
											      matrix=blosum_62, gap_penalty=0, norm=self_score))
		else:
			anchor_1_diff_list = []
			tilted_sequences_list = []
			tilted_template_sequences_list = []
			for i, template_sequence in enumerate(template_sequences):
				alignments = pairwise2.align.localds(self.sequence, template_sequence, blosum_62, -1000, -5)
				if len(alignments) == 0:
					temp_sequence_in_question = ''
					temp_template_sequence = ''
					score_list.append(-1000)
					(temp_anchor_1, temp_anchor_2) = (2, 9)
				else: # When multiple results, return multiple alignment scores.
					max_score = -float("inf")
					best_sequence = ''
					best_template_sequence = ''
					for alignment in alignments:
						temp_sequence_in_question = alignment.seqA
						temp_template_sequence = alignment.seqB
						temp_score = score_sequences(temp_sequence_in_question, temp_template_sequence,
													 anchor_1_list[i], anchor_2_list[i], 
											         matrix=blosum_62, gap_penalty=-0.26, norm=self_score)
						if max_score < temp_score:
							max_score = temp_score
							best_sequence = temp_sequence_in_question
							best_template_sequence = temp_template_sequence
					temp_sequence_in_question = best_sequence
					temp_template_sequence = best_template_sequence
					score_list.append(max_score)
					(temp_anchor_1, temp_anchor_2) = calculate_anchors_given_alignment(temp_sequence_in_question, temp_template_sequence, anchor_1_list[i], anchor_2_list[i])
				anchor_1_diff_list.append(anchor_1_list[i] - temp_anchor_1)
				tilted_sequences_list.append(temp_sequence_in_question)
				tilted_template_sequences_list.append(temp_template_sequence)
			templates['Anchor_diff_1'] = anchor_1_diff_list
			templates['Anchor_diff_2'] = [0] * len(anchor_1_diff_list)
			templates['Tilted_sequence'] = tilted_sequences_list
			templates['Tilted_template_sequence'] = tilted_template_sequences_list

		templates['Peptide_similarity'] = score_list

		templates['Similarity_score'] = 0.5*templates['MHC_similarity'] + 0.5*templates['Peptide_similarity']
		templates = templates.sort_values(by=['Similarity_score'], ascending=False).dropna()

		# 2c. Try filtering by organism first; if that does not bring options, try all!	
		if not templates[templates['MHC'].str[:3] == receptor_allotype[:3]].empty:
			templates_filt = deepcopy(templates[templates['MHC'].str[:3] == receptor_allotype[:3]])
		else:
			templates_filt = deepcopy(templates)
			
		# 2d. Overall similarity (0.5 is empirical, might make it a parameter, we'll see...)
		final_selection = templates_filt[templates_filt['Similarity_score'] > 0]
		cumulative_sum = final_selection['Similarity_score'].sum()
		print("CUMMM: ", cumulative_sum)
		final_selection = final_selection[final_selection['Similarity_score'].cumsum() <= similarity_threshold]
		if (final_selection.empty) or (cumulative_sum < similarity_threshold):
			final_selection = templates[templates['Similarity_score'] > 0] 
			final_selection = final_selection[final_selection['Similarity_score'].cumsum() <= similarity_threshold]
			if final_selection.empty:
				print("No appropriate template an be found for the given peptide-MHC input! Aborting...")
				sys.exit(0)	

		print(final_selection[['pdb_code', 'peptide', 'MHC', 'MHC_similarity', 'Peptide_similarity', 'Similarity_score']])

		return final_selection, anchors, anchor_status

	def initialize_peptide_template(self, template_entry, anchors, anchor_status):

		int_anchors = [int(pos) for pos in anchors.split(",")]
		sequence_length = len(self.sequence)

		peptide_template_file = './new_templates_final/' + template_entry['pdb_code'].values[0]
		template_peptide_length = template_entry['peptide_length'].values[0]

		# 3) Calculate anchor alignment, give the anchor differences that were found previously
		if anchor_status == "Known":
			tilted_sequence, template_tilted_sequence = anchor_alignment(self.sequence, 
																		 template_entry['peptide'].values[0], 
															         	 template_entry['Anchor_diff_1'].values[0], 
															         	 template_entry['Anchor_diff_2'].values[0])
		else:
			tilted_sequence, template_tilted_sequence = template_entry['Tilted_sequence'].values[0], template_entry['Tilted_template_sequence'].values[0]

		# 4) Extract the anchor position numbers for the anchor tolerance step!
		# CAUTION: This could cause inconsistences if the peptide sizes differ greatly, but not really, just making the anchor tolerance step a little bit more obscure
		template_major_anchors = [int(anchor) for anchor in template_entry['Major_anchors'].apply(lambda x: x.split(",")).values[0]]
		template_secondary_anchors = [int(anchor) for anchor in template_entry['Secondary_anchors'].apply(lambda x: x.split(",")).values[0]]
		if anchor_status == "Known": 
			peptide_primary_anchors = int_anchors
		else:
			peptide_primary_anchors = template_major_anchors
			peptide_primary_anchors[1] = template_major_anchors[1] - (len(tilted_sequence) - len(tilted_sequence.rstrip('-'))) # This to adjust the C-terminus anchor when the template is larger in size
			peptide_primary_anchors = [max(1, anchor - template_entry['Anchor_diff_1'].values[0]) for anchor in template_major_anchors] # Anchor adjustement step (see below)

		# 5) Anchors adjustment!
		# Filtering Anchors that won't make sense give the left/right tilt
		Anchor_diff_1 = template_entry['Anchor_diff_1'].values[0]
		peptide_second_anchors = [anchor - Anchor_diff_1 for anchor in template_secondary_anchors]
		peptide_second_anchors = [anchor for anchor in peptide_second_anchors if anchor > 0 and anchor <= sequence_length]
		template_secondary_anchors = [anchor for anchor in template_secondary_anchors if anchor - Anchor_diff_1 > 0 and int(anchor) - Anchor_diff_1 <= sequence_length]

		# 6) Define the peptide template object
		peptide_template = pMHC(pdb_filename=peptide_template_file, 
								peptide=Peptide.frompdb(pdb_filename=peptide_template_file, 
														primary_anchors=template_major_anchors,
														secondary_anchors=template_secondary_anchors,
														tilted_sequence=template_tilted_sequence))

		# 6) Return both the peptide object, as well as the peptide template that was chosen
		self.pdb_filename = None
		self.primary_anchors = peptide_primary_anchors 
		self.secondary_anchors = peptide_second_anchors
		self.tilted_sequence = tilted_sequence

		return peptide_template

	def perform_PTM(self, filestore):

		# Unfortunately, I have to revert to stupid system calls here, because I cannot call pytms from python
		# Maybe one day...
		log_file = filestore + '/03_PTMed_peptides/PTM.log'
		self.pdb_filename = filestore + "/03_PTMed_peptides/PTMed_" + str(self.index) + ".pdb"
		ptm_indexes = []
		for ptm in self.PTM_list:
			PTM, selection = ptm.split(' ', 1)
			ptm_indexes.append(int(selection))
			call(["pymol -qc ./pymol_scripts/" + PTM + ".pml -- " + self.pdb_filename + " " + selection + " " + self.pdb_filename + " > " + log_file + " 2>&1"], shell=True)

		# Further checking the PTMed pdb to see if it actually is canonical and pymol failed; No way to see this from the .log file...
		# When we have no PTM, residue list should be empty and pymol_failed always False.
		ppdb = PandasPdb()
		ppdb.read_pdb(self.pdb_filename)
		pdb_df = ppdb.df['ATOM']
		residues = list(pd.unique(pdb_df[pdb_df['residue_number'].isin(ptm_indexes)]['residue_name']))
		pymol_failed = any(residue in list(standard_three_to_one_letter_code.keys()) for residue in residues)
		if pymol_failed == True:
			with open(filestore + "/05_per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(self.index) + ",Pymol PTM failed!,-\n")
			return True	

		# A. Chain A. I want to make it into chains C as before
		apply_function_to_file(replace_chains, self.pdb_filename, chain_from="A", chain_to="C")

		# B. Weird H0 pymol hydrogens that I want to delete. This is added to other residues during the PTM, so I need to remove them
		apply_function_to_file(delete_elements, self.pdb_filename, element_set=["H0"], chains=["C"])

		# C. I need to re-organize atom indexes, which are a proper mess
		PTMed_tidied = filestore + "/03_PTMed_peptides/PTMed_" + str(self.index) + "tidied.pdb"
		merge_and_tidy_pdb([self.pdb_filename], PTMed_tidied)
		copy_file(PTMed_tidied, self.pdb_filename)
		remove_file(PTMed_tidied)

		return False

	def prepare_for_scoring(self, filestore):
		prep_peptide_loc = "/conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
		self.pdbqt_filename = filestore + "/04_pdbqt_peptides/peptide_" + str(self.index) + ".pdbqt"
		clean = "lps"
		call(["python2.7 " + prep_peptide_loc + " -l " + self.pdb_filename + " -o " + self.pdbqt_filename + " -A None -Z -U " + clean + " -g -s > " + filestore + "/04_pdbqt_peptides/prepare_ligand4.log 2>&1"], shell=True)

		# If the resulting .pdbqt is faulty, delete it. If it does not exist, it is also faulty, so skip whatever else. 
		try:
			seq = pdb_tofasta.run(open(self.pdbqt_filename, 'r'), multi=False)
		except FileNotFoundError:
			with open(filestore + "/05_per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(self.index) + ",Rejected by prepare_ligand4.py,-\n")
			return True

		number_1 = count_number_of_atoms(self.pdb_filename)
		number_2 = count_number_of_atoms(self.pdbqt_filename)
		if(number_1 != number_2):
			remove_file(self.pdbqt_filename)
			with open(filestore + "/05_per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(self.index) + ",Rejected by prepare_ligand4.py,-\n")
			return True
		else:
			return False

	def dock_score_with_SMINA(self, filestore, receptor):

		# SMINA docking and scoring
		self.pdb_filename =  filestore + "/06_scoring_results/model_" + str(self.index) + ".pdb"
		if not receptor.useSMINA and receptor.doMinimization:
			call(["smina -q --scoring vinardo --out_flex " + filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 4 --local_only --minimize --flexres " + receptor.flexible_residues + \
				  " --energy_range 100 --out " + self.pdb_filename + " > " + \
				  filestore + "/06_scoring_results/smina.log 2>&1"], shell=True)
		elif not receptor.useSMINA and not receptor.doMinimization:
			call(["smina -q --scoring vinardo --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 4 --local_only --minimize --energy_range 100 --out " + self.pdb_filename + " > " + \
				  filestore + "/06_scoring_results/smina.log 2>&1"], shell=True)
		elif receptor.useSMINA and receptor.doMinimization:
			call(["smina -q --out_flex " + filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 4 --local_only --minimize --flexres " + receptor.flexible_residues + \
				  " --energy_range 100 --out " + self.pdb_filename + " > " + \
				  filestore + "/06_scoring_results/smina.log 2>&1"], shell=True)
		elif receptor.useSMINA and not receptor.doMinimization:
			call(["smina -q --ligand " + self.pdbqt_filename + \
				  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
				  " --autobox_add 4 --local_only --minimize --energy_range 100 --out " + self.pdb_filename + " > " + \
				  filestore + "/06_scoring_results/smina.log 2>&1"], shell=True)
		atom_number = count_number_of_atoms(self.pdb_filename)
		if(atom_number == 0):
			remove_file(self.pdb_filename)
			with open(filestore + "/05_per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(self.index) + ",Rejected by SMINA (possibly overlapping atoms),-\n")
			return True
		else:
			return False

	def score_with_SMINA(self, filestore, receptor):
		self.pdb_filename = filestore + "/06_scoring_results/model_" + str(self.index) + ".pdb" 
		call(["smina -q --score_only --ligand " + self.pdbqt_filename + \
			  " --receptor " + receptor.pdbqt_filename + " --out " + self.pdb_filename + \
			  " > " + filestore + "/06_scoring_results/smina.log 2>&1"], shell=True)
		move_file(receptor.pdb_filename, filestore + "/09_minimized_receptors/receptor_" + str(self.index) + ".pdb")    

	def compute_anchor_tolerance(self, filestore, receptor, peptide_template_anchors_xyz, anchor_tol, rcd_num_loops):

		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(self.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only anchors
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(self.secondary_anchors)]

		# Only carbon-alpha atoms
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'] == 'CA']

		# Only positions
		pdb_peptide_anchors_xyz = pdb_df_peptide[['x_coord', 'y_coord', 'z_coord']].to_numpy()

		# If difference is smaller than the tolerance, keep the file, else don't (exception for the native one, but it'll probably pass anyway)
		anchor_difference = np.linalg.norm(pdb_peptide_anchors_xyz - peptide_template_anchors_xyz, axis=1)
		if (np.all(anchor_difference < anchor_tol)):
			
			# Keep the scores of the remaining survivors
			with open(self.pdb_filename, 'r') as peptide_handler:
				next(peptide_handler) # Skipping first line
				affinity = peptide_handler.readline().replace("\n", "").split(" ")[2]
			with open(filestore + "/05_per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(self.index) + ",Successfully Modeled," + str(affinity) + "\n")

			dst = filestore + "/08_anchor_filtering/peptide_" + str(self.index) + ".pdb"
			copy_file(self.pdb_filename, dst)
			self.pdb_filename = dst
			return False

		else:
			# delete the minimized receptor coming from SMINA
			if(receptor.doMinimization and os.path.exists(filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb")): remove_file(filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb")
			
			# Keep a log for the anchor difference
			with open(filestore + "/08_anchor_filtering/peptide_" + str(self.index) + ".log", 'a+') as anchor_log:
				anchor_log.write(str(self.index) + "," + ','.join(map(str, anchor_difference))) 
			
			# Keep this result for final printing
			faulty_positions = (anchor_difference > anchor_tol)*self.secondary_anchors
			faulty_positions = " and ".join(np.char.mod('%d', faulty_positions[faulty_positions != 0]))
			with open(filestore + "/05_per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(self.index) + ",Anchor tolerance violated in positions " + faulty_positions + ",-\n")
			
			return True

	def fix_flexible_residues(self, filestore, receptor):

		# Make the flexible receptor output from the SMINA --out_flex argument
		#minimized_receptor_loc = filestore + "/4_SMINA_data/09_minimized_receptors/receptor_" + str(self.index) + ".pdb"
		#if receptor.doMinimization:
		#	call(["python ./helper_scripts/makeflex.py " + \
		#		  filestore + "/4_SMINA_data/receptor_for_smina.pdb " + \
		#		  filestore + "/4_SMINA_data/07_flexible_receptors/receptor_" + str(self.index) + ".pdb " + \
		#		  minimized_receptor_loc],
		#		  stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'), shell=True)

		# Alternative scenario as makeflex.py is probably unstable: Solve the CSP using the CONECT fields to determine the true identity of the atoms

		# Making the CONECT list first:
		edge_list = extract_CONECT_from_pdb(filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb")

		original_ppdb = PandasPdb()
		original_ppdb.read_pdb(filestore + "/receptor_for_smina.pdb")
		original_pdb_df = original_ppdb.df['ATOM']

		flexible_ppdb = PandasPdb()
		flexible_ppdb.read_pdb(filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb")
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

			# print(CA_loc, C_loc)
			matching = csp_solver(sub_edge_list, residue, atom_indexes, CA_loc, C_loc)
			# print(matching)
			# input()
			if matching.shape[0] == 0: # Empty Solution
				# A solution was not found: Most probable case is that the CONECT fields are also broken, meaning that the conformation is invalid as it is. 
				# os.remove(filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb") ## DON'T DELETE FOR KNOW, IN CASE WE HAVE THIS ISSUE AGAIN, INSPECT THE OUTPUT
				with open(filestore + "/05_per_peptide_results/peptide_" + str(self.index) + ".log", 'w') as flexible_log:
					flexible_log.write(str(self.index) + ",Flexible receptor conformation received was faulty,-\n") 
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
		flexible_ppdb.to_pdb(path=filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb", records=['ATOM'], gz=False, append_newline=True)
		minimized_receptor_loc = filestore + "/09_minimized_receptors/receptor_" + str(self.index) + ".pdb"
		merge_and_tidy_pdb([filestore + "/temp_" + str(self.index) + ".pdb", 
							filestore + "/07_flexible_receptors/receptor_" + str(self.index) + ".pdb"],
							minimized_receptor_loc)
		remove_file(filestore + "/temp_" + str(self.index) + ".pdb")

		return False

	def create_peptide_receptor_complexes(self, filestore, receptor):

		# Unify peptide and receptor together
		pMHC_complex = filestore + "/10_pMHC_complexes/pMHC_" + str(self.index) + ".pdb"
		removed = remove_remarks_and_others_from_pdb(self.pdb_filename)
		overwritten = ''.join(removed)
		with open(self.pdb_filename, 'w') as peptide_handler:
			peptide_handler.write(overwritten)
		
		if not receptor.doMinimization: copy_file(filestore + "/receptor_for_smina.pdb",
												  filestore + "/09_minimized_receptors/receptor_" + str(self.index) + ".pdb")
		
		merge_and_tidy_pdb([filestore + "/09_minimized_receptors/receptor_" + str(self.index) + ".pdb", 
							self.pdb_filename], pMHC_complex)


def AA_error_checking(amino_acid):
	if (amino_acid not in standard_three_to_one_letter_code.values()) and (amino_acid not in non_standard_three_to_one_letter_code.values()):
		if verbose(): print("The provided amino acid in the sequence is wrong")
		sys.exit(0)

def process_anchors(anchors, pep_seq):
	# Returns the set of anchors
	pep_length = len(re.sub('[a-z]', '', pep_seq)) # Remove any PTMs that may still exist in the sequence
	anchor_not = sorted(list([anchor_dictionary[str(pep_length)][str(aa_index)] for aa_index in anchors.split(",")]))
	return anchor_not

## PTMs

# Different PTMs

phosphorylation_list = ['pS', 'pT', 'pY']
acetylation_list = ['aK'] # Check details on pytms -> This is not working well, it renames all hydrogens to PyMOL ones
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
			if verbose(): print("PTM chosen is cysteine oxidation. You can only cysteine-hydroxidate (xh), cysteine-oxidate (xo) or cysteine-dioxidate (xd).")
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
