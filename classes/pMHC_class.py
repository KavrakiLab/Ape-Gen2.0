import pymol2

from helper_scripts.Ape_gen_macros import initialize_dir, move_batch_of_files, merge_and_tidy_pdb, all_one_to_three_letter_codes

from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np

from pdbtools import pdb_splitmodel

from subprocess import call
import shutil
import os

from openmm.app import PDBFile, ForceField, Modeller, CutoffNonPeriodic

# PDBFIXER
from pdbfixer import PDBFixer
from openmm.app import PDBFile

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

	def add_sidechains(self, filestore):
		fixer = PDBFixer(filename=self.pdb_filename)
		fixer.findMissingResidues()
		fixer.removeHeterogens(True) #  True keeps water molecules while removing all other heterogens, REVISIT!
		fixer.findMissingAtoms()
		fixer.addMissingAtoms()
		PDBFile.writeFile(fixer.topology, fixer.positions, open(self.pdb_filename, 'w'))

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

		# RCD config -> I think I must include this for residue replacement to work and for no other reason
		# These are also the atoms that I am playing with in RCD (I don't need any others)
		pdb_df_peptide = pdb_df_peptide[pdb_df_peptide['atom_name'].isin(['N', 'CA', 'C', 'O', 'CB'])]

		# Here, I am replacing the residues of the template peptide with the residues of the given peptide.
		for i in range(0, len(pep_seq)):
			pdb_df_peptide.loc[pdb_df_peptide['residue_number'] == i + 1, 'residue_name'] = all_one_to_three_letter_codes[pep_seq[i]]
		
		# Store the peptide now:
		ppdb_peptide.df['ATOM'] = pdb_df_peptide
		anchor_pdb = filestore + '/input_to_RCD/peptide.pdb'
		ppdb_peptide.to_pdb(path=anchor_pdb, records=['ATOM'], gz=False, append_newline=True)

		# Finally, merge those two to create the anchored MHC (peptide contains only the anchors)
		anchored_MHC_file_name = filestore + '/input_to_RCD/anchored_pMHC.pdb'
		merge_and_tidy_pdb([self.pdb_filename, anchor_pdb], anchored_MHC_file_name)
		self.pdb_filename = anchored_MHC_file_name

		# We also have to store the N-terminus and the C-terminus of the peptide for the refinement
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

	def minimizeConf(self, filestore, best_energy, device='CPU'):

		# Read PDB
		pdb = PDBFile(self.pdb_filename)
		top = pdb.getTopology()
		positions = np.array(pdb.positions) 
		numAtoms = len(positions)
		positions = np.reshape(positions, (3*numAtoms,1))

		# Create the ForceField
		forcefield = ForceField('amber/ff14SB.xml', 'amber/phosaa14SB.xml')
		modeller = Modeller(pdb.topology, pdb.positions)
		system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=None)

		# Adding Forces?
		force_constant = 5000
		force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
		force.addGlobalParameter("k", force_constant)
		force.addPerParticleParameter("x0")
		force.addPerParticleParameter("y0")
		force.addPerParticleParameter("z0")
		protein_particles = md.load(filename).top.select("backbone")
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
			r = PDBReporter(filestore + '/OpenMM_confs/' + file, 1)
			r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True))
   		
		return best_energy 