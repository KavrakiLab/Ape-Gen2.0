import pymol2

from helper_scripts.Ape_gen_macros import apply_function_to_file, remove_file, initialize_dir,	   \
											move_batch_of_files, merge_and_tidy_pdb,			   \
											all_one_to_three_letter_codes, replace_CONECT_fields,  \
											merge_connect_fields, select_models, move_file, remove_file, verbose

from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np

from pdbtools import pdb_splitmodel

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

	def align(self, reference, filestore):
		initialize_dir(filestore + '/1_alignment_files')
		#pymol.pymol_argv = ['pymol','-c']
		#pymol.finish_launching()
		p1 = pymol2.PyMOL()
		p1.start()
		p1.cmd.load(self.pdb_filename, "mobile")
		p1.cmd.load(reference.pdb_filename, "ref")

		p1.cmd.align("mobile & chain A", "ref & chain A")

		self.pdb_filename = filestore + '/1_alignment_files/receptor.pdb'
		p1.cmd.save(self.pdb_filename, "mobile")
		p1.cmd.save(filestore + '/1_alignment_files/peptide.pdb', "ref")

		# Also store receptor without peptide and keep that on the receptor part
		# CAUTION: If receptor template is a result of homology modelling, the peptide is ignored there, so
		# the receptor will already be without a peptide to begin with. This does not affect this step at all
		# however
		p1.cmd.create("mobile_sans_peptide", "mobile & chain A")
		self.receptor.pdb_filename = filestore + '/1_alignment_files/receptor_sans_peptide.pdb'
		p1.cmd.save(self.receptor.pdb_filename, "mobile_sans_peptide")

		#pymol.cmd.quit()
		p1.stop()


	def prepare_for_RCD(self, reference, peptide, filestore):

		# Function that prepares files for performing RCD:
		# 1. It removes the peptide from the receptor template
		# 2. It extracts the peptide anchors from the peptime template
		# 3. It replaces those anchors with the amino acids of the peptide that we want to model

		initialize_dir(filestore + '/2_input_to_RCD')

		# 1. Delete the peptide from the receptor template:
		ppdb_receptor = PandasPdb()
		ppdb_receptor.read_pdb(self.pdb_filename)
		pdb_df_receptor = ppdb_receptor.df['ATOM']
		ppdb_receptor.df['ATOM'] = ppdb_receptor.df['ATOM'][ppdb_receptor.df['ATOM']['chain_id'] != 'C']
		self.pdb_filename = filestore + '/2_input_to_RCD/receptor.pdb'
		ppdb_receptor.to_pdb(path=self.pdb_filename, records=['ATOM'], gz=False, append_newline=True)

		# 2. Secondly, keep the anchors and the backbone from the peptide pdb
		ppdb_peptide = PandasPdb()
		ppdb_peptide.read_pdb(reference.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only peptide
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['chain_id'] == 'C']
		
		# RCD config -> I think I must include this for residue replacement to work and for no other reason
		# These are also the atoms that I am playing with in RCD (I don't need any others)
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'].isin(['N', 'CA', 'C', 'O', 'CB'])]

		# Here, I am replacing the (anchor) residues of the template peptide with the residues of the given peptide.
		# Note to self: I don't think I need to replace for the other residues, as this is something RCD takes care of
		anchor_1 = reference.peptide.primary_anchors[0]
		if anchor_1 == 1: anchor_1 = 2
		anchor_2 = reference.peptide.primary_anchors[1]
		template_peptide_len = len(reference.peptide.sequence)
		if anchor_2 == template_peptide_len - 1: anchor_2 = template_peptide_len
		for res in range(1, anchor_1 + 1):
			pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == res, 'residue_name'] = all_one_to_three_letter_codes[peptide.sequence[res - 1]]
		for res in range(anchor_2 - template_peptide_len - 1, 1):
			pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == template_peptide_len + res, 'residue_name'] = all_one_to_three_letter_codes[peptide.sequence[len(peptide.sequence) + res - 1]]
		
		# Removing all middle amino-acids (they will be filled by RCD):
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(list(range(1, anchor_1 + 1)) + list(range(anchor_2 - 1, template_peptide_len + 1)))]

		# But need to re-index those with the peptide length that I want to model in order for RCD to work:
		# In order to not mess up re-indexing, I need to keep as a reference the atom indexes, which are unique
		atom_indexes = {}
		for res in range(anchor_2 - template_peptide_len - 1, 1):
			atom_indexes[res] = pdb_df_peptide[pdb_df_peptide['residue_number'] == template_peptide_len + res]['atom_number'].values
		for res in range(anchor_2 - template_peptide_len - 1, 1):	
			pdb_df_peptide.loc[pdb_df_peptide['atom_number'].isin(atom_indexes[res]), 'residue_number'] = len(peptide.sequence) + res

		# Filter out CBs when we have a Glycine (they shouldn't be there):
		pdb_df_peptide = pdb_df_peptide[~((pdb_df_peptide['residue_name'] == 'GLY') & (pdb_df_peptide['atom_name'] == 'CB'))]

		# Store the peptide now:
		ppdb_peptide.df['ATOM'] = pdb_df_peptide
		anchor_pdb = filestore + '/2_input_to_RCD/peptide.pdb'
		ppdb_peptide.to_pdb(path=anchor_pdb, records=['ATOM'], gz=False, append_newline=True)

		# Finally, merge those two to create the anchored MHC (peptide contains only the anchors)
		# We need to rename B to C, because for some reason C becomes B
		anchored_MHC_file_name = filestore + '/2_input_to_RCD/anchored_pMHC.pdb'
		merge_and_tidy_pdb([self.pdb_filename, anchor_pdb], anchored_MHC_file_name)
		self.pdb_filename = anchored_MHC_file_name

		# We also have to store the N-terminus and the C-terminus of the peptide for the refinement
		anchor_1 = peptide.primary_anchors[0]
		if anchor_1 == 1: anchor_1 = 2
		anchor_2 = peptide.primary_anchors[1]
		if anchor_2 == len(peptide.sequence) - 1: anchor_2 = len(peptide.sequence)
		N_terminus = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(list(range(1, anchor_1)))]
		ppdb_peptide.df['ATOM'] = N_terminus
		ppdb_peptide.to_pdb(path=filestore + '/2_input_to_RCD/N_ter.pdb', 
							records=['ATOM'], gz=False, append_newline=True)
		
		C_terminus = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(list(range(anchor_2, len(peptide.sequence) + 1)))]
		ppdb_peptide.df['ATOM'] = C_terminus
		ppdb_peptide.to_pdb(path=filestore + '/2_input_to_RCD/C_ter.pdb', 
							records=['ATOM'], gz=False, append_newline=True)
		# DONE!

	def RCD(self, peptide, RCD_dist_tol, rcd_num_loops, num_loops, loop_score, filestore):
		
		initialize_dir(filestore + '/3_RCD_data')
		pwd = os.getcwd()

		# Create loops.txt file
		os.chdir(filestore + "/2_input_to_RCD/")
		one_end = peptide.primary_anchors[0] + 1
		if one_end == 2: one_end = 3
		other_end = peptide.primary_anchors[1] - 2
		if other_end == len(peptide.sequence) - 3: other_end = other_end + 1
		with open("./loops.txt", 'w') as loops:
			loops.write("anchored_pMHC.pdb " + str(one_end) + " " + str(other_end) + " C " + peptide.sequence[(one_end - 1):(other_end)])
		loops.close()

		# Call RCD
		call(["rcd -e 1 -x " + pwd + "/RCD_required_files/dunbrack.bin --energy_file " + pwd + "/RCD_required_files/loco.score -o . -d " + str(RCD_dist_tol) + " -n " + str(rcd_num_loops) + " loops.txt >> ../3_RCD_data/rcd.log 2>&1"], shell=True)

		# Move files to back to destination folder
		move_batch_of_files('./', '../3_RCD_data/', query="anchored_pMHC_")
		move_batch_of_files('./', '../3_RCD_data/', query="results")
		os.chdir("../3_RCD_data/")

		# Score loops if the user wants it to
		if verbose(): print("RCD done!")
		# if loop_score != 'none':
		# 	if verbose(): print("Scoring loops with " + loop_score + "...")
		# 	loop_score_to_korpe = {'KORP' : '5', 'ICOSA' : '1'}
		# 	loop_score_to_korpe = loop_score_to_korpe[loop_score]
		# 	file_for_korpe = {'KORP' : 'korp6Dv1.bin', 'ICOSA' : 'loco.score'}
		# 	file_for_korpe = file_for_korpe[loop_score]
		# 	call([pwd + "/RCD_required_files/korpe ../2_input_to_RCD/anchored_pMHC.pdb --loops anchored_pMHC_closed.pdb --score_file " + pwd + "/RCD_required_files/" + file_for_korpe + " -e " + loop_score_to_korpe + " -o korp_res >> korp.log 2>&1 && awk '{$1=$1};1' korp_res.txt > korp_tmp.txt && mv korp_tmp.txt korp_res.txt"], shell=True)
		# 	korp_res = pd.read_csv("korp_res.txt", sep = " ")
		# 	best_indexes = korp_res[korp_res['#loop'] != 0].sort_values(by=['Energy'])['#loop'].head(num_loops).tolist()
		# 	select_models("./anchored_pMHC_closed.pdb", best_indexes , "./anchored_pMHC_closed_temp.pdb")
		# 	move_file("./anchored_pMHC_closed_temp.pdb", "./anchored_pMHC_closed.pdb")
		# else:
		# 	best_indexes = list(range(1, rcd_num_loops + 1))

		# Split the output into files, as the output .pdb has many models			   
		splitted = pdb_splitmodel.run(pdb_splitmodel.check_input(["./anchored_pMHC_closed.pdb"]), outname="model")
		initialize_dir('./splits')
		move_batch_of_files('./', './splits', query="model")
		os.chdir(pwd)
		
		# return best_indexes

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

	def minimizeConf(self, filestore, best_energy, device='CPU'):

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
		force_constant = 5000
		force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
		force.addGlobalParameter("k", force_constant)
		force.addPerParticleParameter("x0")
		force.addPerParticleParameter("y0")
		force.addPerParticleParameter("z0")
		
		#protein_particles = md.load(filename).top.select("backbone")

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
			r = PDBReporter(filestore + '/5_openMM_conformations/14_minimized_complexes/' + file, 1)
			r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True))
			
		return best_energy 