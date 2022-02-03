import pandas as pd
from biopandas.pdb import PandasPdb
import numpy as np

import sys
import os
import re

from helper_scripts.Ape_gen_macros import all_three_to_one_letter_codes, move_batch_of_files, move_file, merge_and_tidy_pdb, replace_chains, remove_remarks_and_others_from_pdb

# PDBFIXER
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from subprocess import call

from pdbtools import pdb_tofasta, pdb_rplchain

class Peptide(object):

	def __init__(self, pdb_filename, sequence):
		self.sequence = sequence
		self.pdb_filename = pdb_filename
		self.pdbqt_filename = None

	@classmethod
	def frompdb(cls, pdb_filename):
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
		return cls(pdb_filename = pdb_filename, sequence = peptide_sequence)
	
	@classmethod
	def fromsequence(cls, peptide_sequence):
		# Initialize peptide from a sequence -> Fetch template!
		templates = pd.read_csv("./template_files/n-mer-templates.csv")
		sequence_length = len(re.sub('[a-z]', '', peptide_sequence)) # Remove PTMs when fetching the template
		peptide_template = templates[templates['Pep_Length'] == sequence_length]['Template'].values[0]
		return cls(pdb_filename = ('./templates/' + peptide_template), sequence = peptide_sequence)

	def add_sidechains(self, filestore, peptide_index):
		fixer = PDBFixer(filename=self.pdb_filename)
		fixer.findMissingResidues()
		fixer.removeHeterogens(True) #  True keeps water molecules while removing all other heterogens, REVISIT!
		fixer.findMissingAtoms()
		fixer.addMissingAtoms()
		#fixer.addMissingHydrogens(7.0) # Ask Mauricio about those
		#fixer.addSolvent(fixer.topology.getUnitCellDimensions()) # Ask Mauricio about those
		self.pdb_filename = filestore + '/SMINA_data/PTMed_peptides/PTMed_' + str(peptide_index) + '.pdb'
		PDBFile.writeFile(fixer.topology, fixer.positions, open(self.pdb_filename, 'w'))

	def perform_PTM(self, filestore, peptide_index, PTM_list):
		# Unfortunately, I have to revert to stupid system calls here, because I cannot call pytms from python
		# Maybe one day...
		log_file = filestore + '/SMINA_data/PTMed_peptides/PTM.log'
		for ptm in PTM_list:
			PTM, selection = ptm.split(' ', 1)
			call(["pymol -qc ./pymol_scripts/" + PTM + ".pml -- " + self.pdb_filename + " " + selection + " " + self.pdb_filename + " > " + log_file + " 2>&1"], shell=True)

		# For some reason, after this step, I get peptide .pdb files with chain A.
		# I want to make it into chains C as before:
		rechained = replace_chains(self.pdb_filename, "A", "C")
		overwritten = ''.join(rechained)
		with open(self.pdb_filename, 'w') as PTMed_file:
			PTMed_file.write(overwritten)

	def prepare_for_scoring(self, filestore, peptide_index, current_round):

		prep_peptide_loc = "/conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
		self.pdbqt_filename = filestore + "/SMINA_data/pdbqt_peptides/peptide_" + str(peptide_index) + ".pdbqt"
		call(["python2.7 " + prep_peptide_loc + " -l " + self.pdb_filename + " -o " + self.pdbqt_filename + " -A bond_hydrogens -Z -U nphs > " + filestore + "/SMINA_data/pdbqt_peptides/prepare_ligand4.log 2>&1"], shell=True)

		# If the resulting .pdbqt is faulty, delete it
		seq = pdb_tofasta.run(open(self.pdbqt_filename, 'r'), multi=False)
		seq = ''.join(seq).split("\n")[1]
		if(len(seq) != len(self.sequence)):
			os.remove(self.pdbqt_filename)
			with open(filestore + "/SMINA_data/per_peptide_results/peptide_" + str(peptide_index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(current_round) + "," + str(peptide_index) + ",Rejected by prepare_ligand4.py,-\n")
			return True
		else:
			return False

	def score_with_SMINA(self, filestore, receptor, peptide_index):

		# SMINA scoring
		self.pdb_filename =  filestore + "/SMINA_data/Scoring_results/model_" + str(peptide_index) + ".pdb"
		if not receptor.useSMINA and receptor.doMinimization:
			call(["smina -q --scoring vinardo --out_flex " + filestore + "/SMINA_data/flexible_receptors/receptor_" + str(peptide_index) + ".pdb --ligand " + self.pdbqt_filename + \
        		  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
        		  " --autobox_add 4 --local_only --minimize --flexres " + receptor.flexible_residues + \
        		  " --energy_range 100 --out " + self.pdb_filename + " > " + \
        		  filestore + "/SMINA_data/Scoring_results/smina.log 2>&1"], shell=True)
		elif not receptor.useSMINA and not receptor.doMinimization:
			call(["smina -q --scoring vinardo --ligand " + self.pdbqt_filename + \
        		  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
        		  " --autobox_add 4 --local_only --minimize --energy_range 100 --out " + self.pdb_filename + " > " + \
        		  filestore + "/SMINA_data/Scoring_results/smina.log 2>&1"], shell=True)
			move_file(receptor.pdbqt_filename, filestore + "/SMINA_data/receptor_smina_min.pdb")
		elif receptor.useSMINA and receptor.doMinimization:
			call(["smina -q --out_flex " + filestore + "/SMINA_data/flexible_receptors/receptor_" + str(peptide_index) + ".pdb --ligand " + self.pdbqt_filename + \
        		  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
        		  " --autobox_add 4 --local_only --minimize --flexres " + receptor.flexible_residues + \
        		  " --energy_range 100 --out " + self.pdb_filename + " > " + \
        		  filestore + "/SMINA_data/Scoring_results/smina.log 2>&1"], shell=True)
		elif receptor.useSMINA and not receptor.doMinimization:
			call(["smina -q --ligand " + self.pdbqt_filename + \
        		  " --receptor " + receptor.pdbqt_filename + " --autobox_ligand " + self.pdbqt_filename + \
        		  " --autobox_add 4 --local_only --minimize --energy_range 100 --out " + self.pdb_filename + " > " + \
        		  filestore + "/SMINA_data/Scoring_results/smina.log 2>&1"], shell=True)
			move_file(receptor.pdbqt_filename, filestore + "/SMINA_data/receptor_smina_min.pdb")		

	def compute_anchor_tolerance(self, filestore, receptor, peptide_template_anchors_xyz, anchor_tol, peptide_index, current_round):

		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(self.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only anchors
		anchors = np.array([1, 2, len(self.sequence) - 1, len(self.sequence)])
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(anchors)]
		
		# Only carbon-alpha atoms
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'] == 'CA']
		
		# Only positions
		pdb_peptide_anchors_xyz = pdb_df_peptide[['x_coord', 'y_coord', 'z_coord']].to_numpy()

		# If difference is smaller than the tolerance, keep the file, else don't
		anchor_difference = np.linalg.norm(pdb_peptide_anchors_xyz - peptide_template_anchors_xyz, axis = 1)
		if np.all(anchor_difference < anchor_tol):
			dst = filestore + "/SMINA_data/Anchor_filtering/peptide_" + str(peptide_index) + ".pdb"
			move_file(self.pdb_filename, dst)
			self.pdb_filename = dst
			return False
		else:
			# delete the minimized receptor coming from SMINA
			if(receptor.doMinimization): os.remove(filestore + "/SMINA_data/flexible_receptors/receptor_" + str(peptide_index) + ".pdb")
			
			# Keep this result for final printing
			faulty_positions = (anchor_difference < anchor_tol)*anchors
			faulty_positions = " and ".join(np.char.mod('%d', faulty_positions[faulty_positions != 0]))
			with open(filestore + "/SMINA_data/per_peptide_results/peptide_" + str(peptide_index) + ".log", 'w') as peptide_handler:
				peptide_handler.write(str(current_round) + "," + str(peptide_index) + ",Anchor tolerance violated in positions " + faulty_positions + ",-\n")
			
			return True

	def create_peptide_receptor_complexes(self, filestore, receptor, peptide_index, current_round):

		# Keep the scores of the remaining survivors
		with open(self.pdb_filename, 'r') as peptide_handler:
			next(peptide_handler) # Skipping first line
			affinity = peptide_handler.readline().replace("\n", "").split(" ")[2]
		with open(filestore + "/SMINA_data/per_peptide_results/peptide_" + str(peptide_index) + ".log", 'w') as peptide_handler:
			peptide_handler.write(str(current_round) + "," + str(peptide_index) + ",Successfully Modeled," + str(affinity) + "\n")

		# Make the flexible receptor output from the SMINA --out_flex argument
		minimized_receptor_loc = filestore + "/SMINA_data/minimized_receptors/receptor_" + str(peptide_index) + ".pdb"
		if receptor.doMinimization:
			call(["python ./helper_scripts/makeflex.py " + \
				  filestore + "/SMINA_data/receptor_for_smina.pdb " + \
				  filestore + "/SMINA_data/flexible_receptors/receptor_" + str(peptide_index) + ".pdb " + \
				  minimized_receptor_loc],
				  stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'), shell=True)

		# Unify peptide and receptor together
		pMHC_complex = filestore + "/Final_conformations/pMHC_" + str(peptide_index) + ".pdb"
		merge_and_tidy_pdb([minimized_receptor_loc, self.pdb_filename], pMHC_complex)
		removed = remove_remarks_and_others_from_pdb(pMHC_complex)
		overwritten = ''.join(removed)
		with open(pMHC_complex, 'w') as pMHC_complex_handler:
			pMHC_complex_handler.write(overwritten)
