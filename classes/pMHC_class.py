import pymol2

from helper_scripts.Ape_gen_macros import apply_function_to_file, remove_file, initialize_dir,	   \
											move_batch_of_files, merge_and_tidy_pdb,			   \
											all_one_to_three_letter_codes, replace_CONECT_fields,  \
											merge_connect_fields, select_models, move_file, 	   \
											remove_file, verbose, add_missing_residues,            \
											apply_mutations, filter_chains, replace_chains

from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np
import random

from pdbtools import pdb_splitmodel, pdb_selmodel

from subprocess import call
import shutil

# OPENMM
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

class pMHC(object):

	def __init__(self, pdb_filename, peptide=None, receptor=None):
		self.pdb_filename = pdb_filename
		self.peptide = peptide # doesn't need PTMs
		self.receptor = receptor

	def remove_peptide(self, filestore):

		p1 = pymol2.PyMOL()
		p1.start()

		p1.cmd.load(self.pdb_filename, "ref")
		p1.cmd.create("sans_peptide", "ref & chain A")
		self.receptor.pdb_filename = filestore + '/receptor_sans_peptide.pdb'
		p1.cmd.save(self.receptor.pdb_filename, "sans_peptide")

		p1.stop()
	
	def align(self, reference, filestore, template_index):
		new_filestore = filestore + '/1_alignment_files/' + str(template_index)

		p1 = pymol2.PyMOL()
		p1.start()

		p1.cmd.load(reference.pdb_filename, "mobile")
		p1.cmd.load(self.pdb_filename, "ref")

		p1.cmd.align("mobile & chain A", "ref & chain A")

		self.pdb_filename = new_filestore + '/receptor.pdb'
		reference.pdb_filename = new_filestore + '/peptide.pdb'
		p1.cmd.save(self.pdb_filename, "ref")
		p1.cmd.save(reference.pdb_filename, "mobile")

		# Also store receptor without peptide and keep that on the receptor part
		# CAUTION: If receptor template is a result of homology modelling, the peptide is ignored there, so
		# the receptor will already be without a peptide to begin with. This does not affect this step at all
		# however
		p1.cmd.create("sans_peptide", "ref & chain A")
		self.receptor.pdb_filename = new_filestore + '/receptor_sans_peptide.pdb'
		p1.cmd.save(self.receptor.pdb_filename, "sans_peptide")

		p1.stop()

	def prepare_for_RCD(self, reference, peptide, filestore, template_index):

		# Let's try a different version of this that does:
		# 1. Deletion
		# 2. Mutation through PDBFixer
		# 3. Insertion through PDBFixer again

		# Extra steps:
	    # 1. Keep copies of the reference file before feeding it to RCD. The purpose of this is to be close to the template and optimize on it, like PANDORA.

		new_filestore = filestore + '/2_input_to_RCD/' + str(template_index)

		# 1. Deletion routine:

		# 1a. Load the peptide template file
		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(reference.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# 1b. Filter for only peptide to be there
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['chain_id'] == 'C']

		# 1c. Calculate anchor differences between reference and peptide in question
		(sequence_in_question, template_sequence) = (peptide.tilted_sequence, reference.peptide.tilted_sequence)

		#print(sequence_in_question)
		#print(template_sequence)

		indexes_for_deletion = [pos + 1 for pos, char in enumerate(sequence_in_question) if char == '-']
		indexes_for_deletion = [pos - (len(template_sequence) - len(template_sequence.lstrip('-'))) if pos > 5 else pos for pos in indexes_for_deletion] # Adjustment for C-terminus indexes when there is a tilt from the front (example: 2FWO).

		pdb_df_peptide = pdb_df_peptide[~pdb_df_peptide['residue_number'].isin(indexes_for_deletion)]
		pdb_df_peptide['residue_number'] = pdb_df_peptide['residue_number'] - (len(sequence_in_question) - len(sequence_in_question.lstrip('-'))) + (len(template_sequence) - len(template_sequence.lstrip('-')))

		# 1e. Store
		ppdb_peptide.df['ATOM'] = pdb_df_peptide
		anchor_pdb = new_filestore + '/2_peptide_del.pdb'
		ppdb_peptide.to_pdb(path=anchor_pdb, records=['ATOM'], gz=False, append_newline=True)
		anchored_MHC_file_name = new_filestore + '/anchored_pMHC_del.pdb'
		merge_and_tidy_pdb([self.receptor.pdb_filename, anchor_pdb], anchored_MHC_file_name)

		# 2. Mutate

		# 2a. Calculate the mutation list
		aa_list_ref = list(template_sequence)
		aa_list_in_question = list(sequence_in_question)
		mutation_list = []

		i = 0
		j = i + 1
		while(i < len(sequence_in_question)):
			if aa_list_in_question[i] == '-':
				j -= 1
			elif aa_list_ref[i] != '-' and aa_list_in_question[i] != '-' and aa_list_ref[i] != aa_list_in_question[i]:
				mutation_list.append(all_one_to_three_letter_codes[aa_list_ref[i]] + "-" + str(j) + "-" + all_one_to_three_letter_codes[aa_list_in_question[i]])
			else:
				pass
			i += 1
			j += 1

		# 2b. Apply mutations using the mutation list
		apply_mutations(anchored_MHC_file_name, new_filestore, mutation_list)

		# 3. Insert

		# 3a. To make the insert, the SEQRES field needs to be added on top of the file
		# See documentation: https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/primary-sequences-and-the-pdb-format
		# Obviously, this is a pain, because depending on the size of the peptide, there are different conventions:
		seqres_pdb = new_filestore + '/seqres.pdb'
		three_letter_list = [all_one_to_three_letter_codes[aa] for aa in list(peptide.sequence)]
		if len(three_letter_list) <= 9:
			three_letter_string = ' '.join(three_letter_list)
			with open(seqres_pdb, 'w') as seqres:
				seqres.write("SEQRES   1 C    " + str(len(three_letter_list)) + "  " + three_letter_string)
			seqres.close()
		elif len(three_letter_list) <= 13:
			three_letter_string = ' '.join(three_letter_list)
			with open(seqres_pdb, 'w') as seqres:
				seqres.write("SEQRES   1 C   " + str(len(three_letter_list)) + "  " + three_letter_string)
			seqres.close()
		else:
			three_letter_string_1 = ' '.join(three_letter_list[:13])
			three_letter_string_2 = ' '.join(three_letter_list[13:])
			with open(seqres_pdb, 'w') as seqres:
				seqres.write("SEQRES   1 C   " + str(len(three_letter_list)) + "  " + three_letter_string_1)
				seqres.write("\nSEQRES   2 C   " + str(len(three_letter_list)) + "  " + three_letter_string_2)
			seqres.close()
		anchored_MHC_file_name_seqres = new_filestore + '/anchored_pMHC_native.pdb'
		merge_and_tidy_pdb([seqres_pdb, anchored_MHC_file_name], anchored_MHC_file_name_seqres)

		# 3b. Inserting the extra amino acids with PDBFixer:
		add_missing_residues(anchored_MHC_file_name_seqres, new_filestore)

		# Keep copies of the peptide file separately, as it will be used for refinement downstream:
		filter_chains(anchored_MHC_file_name_seqres, ("C",), new_filestore + "/model_" + str(template_index) + ".pdb")

		# Keep only the backbone in the end:
		ppdb_peptide.read_pdb(anchored_MHC_file_name_seqres)
		pdb_df_peptide = ppdb_peptide.df['ATOM'].copy()
		pdb_df_peptide = pdb_df_peptide[(pdb_df_peptide['chain_id'] == 'A') | \
										((pdb_df_peptide['chain_id'] == 'C') & \
										(pdb_df_peptide['atom_name'].isin(['N', 'CA', 'C', 'O', 'CB'])))]

		# Filter out CBs when we have a Glycine (they shouldn't be there):
		pdb_df_peptide = pdb_df_peptide[~((pdb_df_peptide['residue_name'] == 'GLY') & \
										 (pdb_df_peptide['atom_name'] == 'CB') & \
										 ((pdb_df_peptide['chain_id'] == 'C')))]
		
		min_residue_number = pdb_df_peptide['residue_number'].min()
		pdb_df_peptide.update(pdb_df_peptide['residue_number'] + 1 - min_residue_number) 

		ppdb_peptide.df['ATOM'] = pdb_df_peptide
		anchored_MHC_file_name = new_filestore + '/anchored_pMHC_proper.pdb'								
		ppdb_peptide.to_pdb(path=anchored_MHC_file_name, records=['ATOM'], gz=False, append_newline=True)

		# We also have to store the N-terminus and the C-terminus of the peptide for the refinement
		anchor_1 = peptide.primary_anchors[0]
		if anchor_1 == 1: anchor_1 = 2
		anchor_2 = peptide.primary_anchors[1]

		# Only peptide
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['chain_id'] == 'C']

		N_terminus = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(list(range(1, anchor_1)))]
		ppdb_peptide.df['ATOM'] = N_terminus
		ppdb_peptide.to_pdb(path=new_filestore + '/N_ter.pdb', 
							records=['ATOM'], gz=False, append_newline=True)
		
		C_terminus = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(list(range(anchor_2, len(peptide.sequence) + 1)))]
		ppdb_peptide.df['ATOM'] = C_terminus
		ppdb_peptide.to_pdb(path=new_filestore + '/C_ter.pdb', 
							records=['ATOM'], gz=False, append_newline=True)
		# DONE!

	def RCD(self, peptide, RCD_dist_tol, rcd_num_loops, filestore, template_index):

		pwd = os.getcwd()
		os.chdir(filestore + "/2_input_to_RCD/" + str(template_index) + "/")
		
		# Create loops.txt file
		# The canonical case is the N-termini anchor being in pos 2, however, if it is in pos 1, just because we need
		# to have a Start-2 Start-1 config, the N-termini endpoint must be in pos 3
		one_end = peptide.primary_anchors[0] + 1
		if one_end == 2: one_end = 3

		# Similarly, the canonical case in C-termini is P-Omega. However, in that case, we need an extra AA because the
		# Configuration is End-2 End-1. So in the canonical case, P-Omega -2 is the endpoint. 
		other_end = peptide.primary_anchors[1] - 2
		with open("./loops.txt", 'w') as loops:
			loops.write("anchored_pMHC_proper.pdb " + str(one_end) + " " + str(other_end) + " C " + peptide.sequence[(one_end - 1):(other_end)])
		loops.close()

		# Call RCD
		call(["rcd -e 1 -x " + pwd + "/RCD_required_files/dunbrack.bin --energy_file " + pwd + "/helper_files/loco.score -o . -d " + str(RCD_dist_tol) + " -n " + str(rcd_num_loops) + " --bench loops.txt >> ../../3_RCD_data/" + str(template_index) + "/rcd.log 2>&1 && awk '{$1=$1};1' anchored_pMHC_proper_rmsd.txt > korp_tmp.txt && mv korp_tmp.txt anchored_pMHC_proper_rmsd.txt"], shell=True)

		# Move files to back to destination folder
		move_batch_of_files('./', '../../3_RCD_data/' + str(template_index) + '/', query="anchored_pMHC_proper_")
		move_batch_of_files('./', '../../3_RCD_data/' + str(template_index) + '/', query="results")
		os.chdir(pwd)

		# DONE!
		
	def loop_ranking(self, rcd_num_loops, num_loops, loop_score, non_sampled_confs, filestore, template_index):
		
		pwd = os.getcwd()
		os.chdir(filestore + '/3_RCD_data/' + str(template_index) + '/')

		korp_res = pd.read_csv("anchored_pMHC_proper_rmsd.txt", sep = " ", comment='#', header=None)
		korp_res.columns = ['Loop', 'RMSD', 'Bump', 'BumpEx', 'BumpIn', 'Energy']
		if loop_score == 'KORP':
			call([pwd + "/RCD_required_files/korpe ../../2_input_to_RCD/" + str(template_index) + "/anchored_pMHC_proper.pdb --loops anchored_pMHC_proper_closed.pdb --score_file " + pwd + "/RCD_required_files/korp6Dv1.bin -e 5 -o korp_res >> korp.log 2>&1 && awk '{$1=$1};1' korp_res.txt > korp_tmp.txt && mv korp_tmp.txt korp_res.txt"], shell=True)
			korp_res = pd.read_csv("korp_res.txt", sep = " ", comment='#', header=None)
			korp_res.columns = ['Loop', 'Energy']
			best_indexes = korp_res[korp_res['Loop'] != 0].sort_values(by=['Energy'])['Loop'].head(num_loops).astype(int).tolist()
		elif loop_score == 'ICOSA':
			best_indexes = korp_res.sort_values(by=['Energy'])['Loop'].head(num_loops).astype(int).tolist()
		elif loop_score == 'RMSD':
			best_indexes = korp_res.sort_values(by=['RMSD'])['Loop'].head(num_loops).astype(int).tolist()
		else:
			best_indexes = random.sample(range(rcd_num_loops), num_loops)
		best_indexes = [idx + 1 for idx in best_indexes]

		# Split the output into files, as the output .pdb has many models			   
		splitted = pdb_splitmodel.run(pdb_splitmodel.check_input(["./anchored_pMHC_proper_closed.pdb"]), outname="model")
		initialize_dir('./splits')
		move_batch_of_files('./', './splits', query="model")
		os.chdir(pwd)
		
		# Add the rest of the non_sampled conformations, same way as PANDORA does. 
		best_indexes.extend(list(range(rcd_num_loops + 1, rcd_num_loops + non_sampled_confs + 1)))

		return best_indexes

	def set_anchor_xyz(self, anchor_selection, peptide):

		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(self.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only peptide
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['chain_id'] == 'C'] 

		# Only anchors (based on selection)
		if anchor_selection == 'primary':
			pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(self.peptide.primary_anchors)]
			tolerance_anchors = peptide.primary_anchors
		elif anchor_selection == 'secondary':
			pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(self.peptide.secondary_anchors)]
			tolerance_anchors = peptide.secondary_anchors
		else:
			pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin([])]
			tolerance_anchors = []

		# Only carbon-alpha atoms
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'] == 'CA']
		
		# Only positions
		pdb_df_peptide = pdb_df_peptide[['x_coord', 'y_coord', 'z_coord']]

		return pdb_df_peptide.to_numpy(), tolerance_anchors

	def add_PTM_CONECT_fields(self, filestore, PTM_list, peptide_index):

		ppdb_complex = PandasPdb()
		ppdb_complex.read_pdb(self.pdb_filename)
		pdb_df_complex = ppdb_complex.df['ATOM']

		# Only peptide
		pdb_df_complex = pdb_df_complex[pdb_df_complex['chain_id'] == 'C'] 

		file_list = [self.pdb_filename]
		for PTM in PTM_list:
			PTM_index = int(PTM.split(' ')[1])
			sub_pdb = pdb_df_complex[pdb_df_complex['residue_number'] == PTM_index]
			residue_name = pd.unique(sub_pdb['residue_name'])[0]

			# Define external bonds:
			external_bonds_list = []
			if PTM_index != 1:
				previous_C = pdb_df_complex[(pdb_df_complex['residue_number'] == PTM_index - 1) & (pdb_df_complex['atom_name'] == 'C')]['atom_number'].item()
				current_N = sub_pdb[sub_pdb['atom_name'] == 'N']['atom_number'].item()
				external_bonds_list.append((previous_C, current_N))
			if PTM_index != len(self.peptide.sequence):
				current_C = sub_pdb[sub_pdb['atom_name'] == 'C']['atom_number'].item()
				next_N = pdb_df_complex[(pdb_df_complex['residue_number'] == PTM_index + 1) & (pdb_df_complex['atom_name'] == 'N')]['atom_number'].item()
				external_bonds_list.append((current_C, next_N))

			conect_file = apply_function_to_file(func=replace_CONECT_fields, 
												 input_filename='./PTM_residue_templates/' + residue_name + '.conect', 
												 output_filename=filestore + '/5_openMM_conformations/12_PTM_conect_indexes/conect_' + 
												 					str(peptide_index)  + residue_name + str(PTM_index) + '.pdb', 
												 index_df=sub_pdb,
												 external_bonds_list=external_bonds_list)
			file_list.append(conect_file)


		self.pdb_filename = filestore + '/5_openMM_conformations/13_connected_pMHC_complexes/pMHC_' + str(peptide_index) + '.pdb'
		merge_connect_fields(file_list, self.pdb_filename)

	def minimizeConf(self, filestore, best_energy, no_constraints_openmm=False, device='CPU'):

		# Read PDB
		pdb = PDBFile(self.pdb_filename)
		top = pdb.getTopology()
		positions = np.array(pdb.positions) 
		numAtoms = len(positions)
		positions = np.reshape(positions, (3*numAtoms,1))

		# Create the ForceField
		forcefield = ForceField('amber/ff14SB.xml', 'amber/phosaa14SB.xml')
		# forcefield = ForceField('charmm/charmm36_nowaters.xml')
		modeller = Modeller(pdb.topology, pdb.positions)
		system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=None)

		# Adding Forces?
		if not no_constraints_openmm:
			force_constant = 5000
			force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
			force.addGlobalParameter("k", force_constant)
			force.addPerParticleParameter("x0")
			force.addPerParticleParameter("y0")
			force.addPerParticleParameter("z0")

			ppdb = PandasPdb()
			ppdb.read_pdb(self.pdb_filename)
			pdb_df = ppdb.df['ATOM']
			protein_particles = [atomind - 1 for atomind in pdb_df[pdb_df['atom_name'].isin(["N", "O", "C", "CA"])]['atom_number'].tolist()]
			particle_indices = []
			for protein_particle in protein_particles:
				particle_indices.append(force.addParticle(int(protein_particle), modeller.positions[protein_particle]) )
			system.addForce(force)

		# Enumerate forces?
		forces = system.getForces()
		for i, f in enumerate(forces):
			f.setForceGroup(i)

		# Rest (Integrator + Platform)
		integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
		platform = Platform.getPlatformByName(device)

		# Create Simulation
		simulation = Simulation(modeller.topology, system, integrator, platform)
		simulation.context.setPositions(modeller.positions)

		# Minimize energy
		simulation.minimizeEnergy()
		simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, potentialEnergy=True, 
									temperature=True, progress=False, remainingTime=True, speed=True, 
									totalSteps=250000, separator='\t'))

		# Write results to a new file if energy is small enough
		energy = simulation.context.getState(getEnergy=True).getPotentialEnergy() / kilojoule_per_mole
		if energy < best_energy:
			best_energy = energy
			path, file = os.path.split(self.pdb_filename)
			r = PDBReporter(filestore + '/5_openMM_conformations/10_pMHC_complexes/' + file, 1)
			r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True))
		
		# Rename the chains, for some reason it's chain B and not C
		apply_function_to_file(replace_chains, filestore + '/5_openMM_conformations/10_pMHC_complexes/' + file, chain_from="B", chain_to="C")	

		return best_energy 
