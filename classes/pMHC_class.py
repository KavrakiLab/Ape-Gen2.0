import pymol2

<<<<<<< HEAD
from helper_scripts.Ape_gen_macros import apply_function_to_file, remove_file, initialize_dir,	   \
											move_batch_of_files, merge_and_tidy_pdb,			   \
											all_one_to_three_letter_codes, replace_CONECT_fields,  \
											merge_connect_fields, verbose
=======
from helper_scripts.Ape_gen_macros import remove_file, initialize_dir, move_batch_of_files, merge_and_tidy_pdb, all_one_to_three_letter_codes, replace_CONECT_fields, merge_connect_fields
>>>>>>> made os.remove a macro

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

# PDBFIXER
from pdbfixer import PDBFixer
from openmm.app import PDBFile

class pMHC(object):

	def __init__(self, pdb_filename, peptide = None, receptor = None):
		self.pdb_filename = pdb_filename
		self.peptide = peptide # doesn't need PTMs
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

		p1.cmd.align("mobile & chain A", "ref & chain A")

		self.pdb_filename = filestore + '/alignment_files/receptor.pdb'
		p1.cmd.save(self.pdb_filename, "mobile")
		p1.cmd.save(filestore + '/alignment_files/peptide.pdb', "ref")

		# Also store receptor without peptide and keep that on the receptor part
		# CAUTION: If receptor template is a result of homology modelling, the peptide is ignored there, so
		# the receptor will already be without a peptide to begin with. This does not affect this step at all
		# however
		p1.cmd.create("mobile_sans_peptide", "mobile & chain A")
		self.receptor.pdb_filename = filestore + '/alignment_files/receptor_sans_peptide.pdb'
		p1.cmd.save(self.receptor.pdb_filename, "mobile_sans_peptide")

		#pymol.cmd.quit()
		p1.stop()

	def add_sidechains(self, filestore):
		fixer = PDBFixer(filename=self.pdb_filename)
		fixer.findMissingResidues()
		fixer.removeHeterogens(True) #  True keeps water molecules while removing all other heterogens, REVISIT!
		fixer.findMissingAtoms()
		fixer.addMissingAtoms()
		PDBFile.writeFile(fixer.topology, fixer.positions, open(self.pdb_filename, 'w'), keepIds=True)

	def prepare_for_RCD(self, reference, filestore):

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
		
		# RCD config -> I think I must include this for residue replacement to work and for no other reason
		# These are also the atoms that I am playing with in RCD (I don't need any others)
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'].isin(['N', 'CA', 'C', 'O', 'CB'])]

		# It is important here to ensure that the peptide template length and the peptide tomodel length are equal
		# If the length of the template is larger, we need to remove amino acids (removing the middle ones)
		# If the length of the template is smaller, we need to add aminoacids

		# Here, I am replacing the (anchor) residues of the template peptide with the residues of the given peptide.
		# Note to self: I don't think I need to replace for the other residues, as this is something RCD takes care of
		template_peptide_len = pdb_df_peptide['residue_number'].max()
		pep_seq = self.peptide.sequence
		pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == 1, 'residue_name'] = all_one_to_three_letter_codes[pep_seq[0]]
		pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == 2, 'residue_name'] = all_one_to_three_letter_codes[pep_seq[1]]
		pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == template_peptide_len - 1, 'residue_name'] = all_one_to_three_letter_codes[pep_seq[len(pep_seq) - 2]]
		pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == template_peptide_len, 'residue_name'] = all_one_to_three_letter_codes[pep_seq[len(pep_seq) - 1]]
		
		# Removing all middle amino-acids (they will be filled by RCD):
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin([1, 2, template_peptide_len - 1, template_peptide_len])]

		# But need to re-index those with the peptide length that I want to model in order for RCD to work:
		# In order to not mess up re-indexing, I need to keep as a reference the atom indexes, which are unique
		atom_indexes_c1 = pdb_df_peptide[pdb_df_peptide['residue_number'] == template_peptide_len - 1]['atom_number'].values
		atom_indexes_c = pdb_df_peptide[pdb_df_peptide['residue_number'] == template_peptide_len]['atom_number'].values
		pdb_df_peptide.loc[pdb_df_peptide['atom_number'].isin(atom_indexes_c1), 'residue_number'] = len(pep_seq) - 1
		pdb_df_peptide.loc[pdb_df_peptide['atom_number'].isin(atom_indexes_c), 'residue_number'] = len(pep_seq)

		# Store the peptide now:
		ppdb_peptide.df['ATOM'] = pdb_df_peptide
		anchor_pdb = filestore + '/input_to_RCD/peptide.pdb'
		ppdb_peptide.to_pdb(path=anchor_pdb, records=['ATOM'], gz=False, append_newline=True)

		# Finally, merge those two to create the anchored MHC (peptide contains only the anchors)
		# We need to rename B to C, because for some reason C becomes B
		anchored_MHC_file_name = filestore + '/input_to_RCD/anchored_pMHC.pdb'
		merge_and_tidy_pdb([self.pdb_filename, anchor_pdb], anchored_MHC_file_name)
		self.pdb_filename = anchored_MHC_file_name

		# We also have to store the N-terminus and the C-terminus of the peptide for the refinement
		N_terminus = pdb_df_peptide[pdb_df_peptide['residue_number'] == 1]
		ppdb_peptide.df['ATOM'] = N_terminus
		ppdb_peptide.to_pdb(path=filestore + '/input_to_RCD/N_ter.pdb', 
							records=['ATOM'], gz=False, append_newline=True)

		C_terminus = pdb_df_peptide[pdb_df_peptide['residue_number'] == len(pep_seq)]
		ppdb_peptide.df['ATOM'] = C_terminus
		ppdb_peptide.to_pdb(path=filestore + '/input_to_RCD/C_ter.pdb', 
							records=['ATOM'], gz=False, append_newline=True)

		# DONE!

	def RCD(self, RCD_dist_tol, num_loops, filestore):

		initialize_dir(filestore + '/RCD_data')

		# Create loops.txt file
		last_non_anchor = len(self.peptide.sequence) - 2
		with open(filestore + "/input_to_RCD/loops.txt", 'w') as loops:
			loops.write(filestore + "/input_to_RCD/anchored_pMHC.pdb 3 " + str(last_non_anchor) + " C " + self.peptide.sequence[2:last_non_anchor])
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
		remove_file(filestore + '/../../results.txt')

	def set_anchor_xyz(self, reference, anchors):

		ppdb_peptide = PandasPdb()

		ppdb_peptide.read_pdb(reference.pdb_filename)
		pdb_df_peptide = ppdb_peptide.df['ATOM']

		# Only peptide
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['chain_id'] == 'C'] 

		# Only anchors
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['residue_number'].isin(anchors)]

		# Only carbon-alpha atoms
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'] == 'CA']
		
		# Only positions
		pdb_df_peptide = pdb_df_peptide[['x_coord', 'y_coord', 'z_coord']]

		return pdb_df_peptide.to_numpy()

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
												 output_filename=filestore + '/OpenMM_confs/PTM_conect_indexes/conect_' + 
												 					str(peptide_index)  + residue_name + str(PTM_index) + '.pdb', 
												 index_df=sub_pdb,
												 external_bonds_list=external_bonds_list)
			file_list.append(conect_file)


		self.pdb_filename = filestore + '/OpenMM_confs/connected_pMHC_complexes/pMHC_' + str(peptide_index) + '.pdb'
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
			r = PDBReporter(filestore + '/OpenMM_confs/minimized_complexes/' + file, 1)
			r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True))
   			
		return best_energy 