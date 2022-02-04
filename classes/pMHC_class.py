import pymol2

from helper_scripts.Ape_gen_macros import initialize_dir, move_batch_of_files, merge_and_tidy_pdb, all_one_to_three_letter_codes

from biopandas.pdb import PandasPdb
import pandas as pd
from pdbtools import pdb_splitmodel

from subprocess import call
import shutil
import os

class pMHC(object):

	def __init__(self, pdb_filename, peptide = None, receptor = None):
		self.pdb_filename = pdb_filename
		self.peptide = peptide
		self.receptor = receptor
		self.anchor_xyz = None

	def align(self, reference, filestore):
		initialize_dir(filestore + '/alignment_files')
		#pymol.pymol_argv = ['pymol','-c']
		#pymol.finish_launching()
		p1 = pymol2.PyMOL()
		p1.start()
		p1.cmd.load(self.pdb_filename, "mobile")
		p1.cmd.load(reference.pdb_filename, "ref")

		p1.cmd.align("mobile & (chain A | chain B)", "ref & (chain A | chain B)")

		self.pdb_filename = filestore + '/alignment_files/receptor.pdb'
		p1.cmd.save(self.pdb_filename, "mobile")
		p1.cmd.save(filestore + '/alignment_files/peptide.pdb', "ref")

		# Also store receptor without peptide and keep that on the receptor part
		p1.cmd.create("mobile_sans_peptide", "mobile & (chain A | chain B)")
		self.receptor.pdb_filename = filestore + '/alignment_files/receptor_sans_peptide.pdb'
		p1.cmd.save(self.receptor.pdb_filename, "mobile_sans_peptide")

		#pymol.cmd.quit()
		p1.stop()

	def prepare_for_RCD(self, reference, filestore, pep_seq):

		# Function that prepares files for performing RCD
		# Namely, it extracts the peptide anchors from the peptime template
		# It removes the peptide from the receptor template
		# It unifies those results, making the receptor + anchors that we want to model using RCD

		initialize_dir(filestore + '/input_to_RCD')

		# First, delete the peptide from the receptor template:
		ppdb_receptor = PandasPdb()
		ppdb_receptor.read_pdb(self.pdb_filename)
		pdb_df_receptor = ppdb_receptor.df['ATOM']
		ppdb_receptor.df['ATOM'] = ppdb_receptor.df['ATOM'][ppdb_receptor.df['ATOM']['chain_id'] != 'C']
		self.pdb_filename = filestore + '/input_to_RCD/receptor.pdb'
		ppdb_receptor.to_pdb(path=self.pdb_filename, records=['ATOM'], gz=False, append_newline=True)

		# Secondly, keep the anchors and the backbone from the peptide pdb
		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(reference.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only peptide
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['chain_id'] == 'C'] 

		# Only anchors
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin([1, 2, len(pep_seq) - 1, len(pep_seq)])]

		# RCD config -> I think I must include this for residue replacement to work and for no other reason
		# These are also the atoms that I am playing with in RCD (I don't need any others)
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'].isin(['N', 'CA', 'C', 'O', 'CB'])]

		# Here, I am replacing the anchor residues of the template peptide with the residues of the given peptide
		pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == 1, 'residue_name'] = all_one_to_three_letter_codes[pep_seq[0]]
		pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == 2, 'residue_name'] = all_one_to_three_letter_codes[pep_seq[1]]
		pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == len(pep_seq) - 1, 'residue_name'] = all_one_to_three_letter_codes[pep_seq[len(pep_seq) - 2]]
		pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == len(pep_seq), 'residue_name'] = all_one_to_three_letter_codes[pep_seq[len(pep_seq) - 1]]

		# Store the anchors now:
		ppdb_peptide.df['ATOM'] = pdb_df_peptide
		anchor_pdb = filestore + '/input_to_RCD/anchors.pdb'
		ppdb_peptide.to_pdb(path=anchor_pdb, records=['ATOM'], gz=False, append_newline=True)

		#Finally, merge those two to create the anchored MHC (peptide contains only the anchors)
		anchored_MHC_file_name = filestore + '/input_to_RCD/anchored_pMHC.pdb'
		merge_and_tidy_pdb([self.pdb_filename, anchor_pdb], anchored_MHC_file_name)
		self.pdb_filename = anchored_MHC_file_name

		# Before we go, we also have to store the N-terminus and the C-terminus of the peptide for the refinement
		N_terminus = pdb_df_peptide[pdb_df_peptide['residue_number'] == 1]
		ppdb_peptide.df['ATOM'] = N_terminus
		ppdb_peptide.to_pdb(path=filestore + '/input_to_RCD/N_ter.pdb', 
							records=['ATOM'], gz=False, append_newline=True)

		C_terminus = pdb_df_peptide[pdb_df_peptide['residue_number'] == len(reference.peptide.sequence)]
		ppdb_peptide.df['ATOM'] = C_terminus
		ppdb_peptide.to_pdb(path=filestore + '/input_to_RCD/C_ter.pdb', 
						    records=['ATOM'], gz=False, append_newline=True)

		# DONE!

	def RCD(self, peptide, RCD_dist_tol, num_loops, filestore):

		initialize_dir(filestore + '/RCD_data')

		# Create loops.txt file
		last_non_anchor = len(peptide.sequence) - 2
		with open(filestore + "/input_to_RCD/loops.txt", 'w') as loops:
			loops.write(filestore + "/input_to_RCD/anchored_pMHC.pdb 3 " + str(last_non_anchor) + " C " + peptide.sequence[2:last_non_anchor])
		loops.close()

		# Perform RCD:
		call(["rcd -e 1 -x ./RCD_required_files/dunbrack.bin --energy_file ./RCD_required_files/loco.score -o . -d " + str(RCD_dist_tol) + " -n " + str(num_loops) + " " + filestore + "/input_to_RCD/loops.txt >> " + filestore + "/RCD_data/rcd.log 2>&1"], shell=True)
 		
 		# Move files to back to destination folder (think about making a function for this)
		move_batch_of_files(filestore + '/input_to_RCD/', filestore + '/RCD_data', query = "anchored_pMHC_")

		# Split the output into files, as the output .pdb has many models		
		splitted = pdb_splitmodel.run(pdb_splitmodel.check_input([filestore + "/RCD_data/anchored_pMHC_closed.pdb"],
																  ), outname = "model")
		initialize_dir(filestore + '/RCD_data/splits')	
		move_batch_of_files('./', filestore + '/RCD_data/splits', query = "model")

	def set_anchor_xyz(self, reference, pep_seq):

		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(reference.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only peptide
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['chain_id'] == 'C'] 

		# Only anchors
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin([1, 2, len(pep_seq) - 1, len(pep_seq)])]

		# Only carbon-alpha atoms
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'] == 'CA']
		
		# Only positions
		pdb_df_peptide = pdb_df_peptide[['x_coord', 'y_coord', 'z_coord']]

		return pdb_df_peptide.to_numpy()
