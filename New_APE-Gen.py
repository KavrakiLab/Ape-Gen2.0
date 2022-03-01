from helper_scripts import argparser
from helper_scripts.Ape_gen_macros import initialize_dir, merge_and_tidy_pdb, sequence_PTM_processing, create_csv_from_list_of_files, copy_file, pretty_print_analytics, move_batch_of_files, copy_batch_of_files

from classes.Peptide_class import Peptide
from classes.Receptor_class import Receptor
from classes.pMHC_class import pMHC

import pandas as pd
import numpy as np
import re

import sys
from mpire import WorkerPool
from tqdm import tqdm

# Temporary
from subprocess import call
from pdbtools import pdb_mkensemble
import glob
from time import sleep
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def peptide_refinement_and_scoring(peptide_index, filestore, PTM_list, receptor, peptide_template_anchors_xyz, anchor_tol, current_round):

	# Routine that refines and scores a peptide/receptor pair with SMINA/Vinardo

	# 1. Assemble peptide by mergin the peptide anchors and the middle part
	model_location = filestore + '/RCD_data/splits/model_' + str(peptide_index) + '.pdb'
	Nterm_location = filestore + '/input_to_RCD/N_ter.pdb'
	Cterm_location = filestore + '/input_to_RCD/C_ter.pdb'
	assembled_peptide = filestore + '/SMINA_data/assembled_peptides/assembled_' + str(peptide_index) + '.pdb'
	merge_and_tidy_pdb([Nterm_location, model_location, Cterm_location], assembled_peptide)
	peptide = Peptide.frompdb(assembled_peptide)

	#peptide_is_not_valid = peptide.prepare_for_scoring(filestore, peptide_index, current_round)
	#if(peptide_is_not_valid): return

	# 2. Now that the peptide is assembled, Fill in the sidechains with pdbfixer
	peptide.add_sidechains(filestore, peptide_index)

	# 3. Do PTMs
	peptide.perform_PTM(filestore, peptide_index, PTM_list)

	# Intermediate OpenMM step that fixes things maybe?
	#peptide.minimizeConf(filestore, peptide_index)
	#input()

	# 4. Score with SMINA
	# 4a. .pdb to .pdbqt transformation using autodocktools routines (very good for filtering bad conformations)
	peptide_is_not_valid = peptide.prepare_for_scoring(filestore, peptide_index, current_round)
	if(peptide_is_not_valid): return

	# 4b. Score with SMINA (or other options, depending on args)
	peptide.score_with_SMINA(filestore, receptor, peptide_index) 

	# 5. Anchor filtering (based on anchor tolerance argument) (improvement from previous version)
	peptide_is_not_valid = peptide.compute_anchor_tolerance(filestore, receptor, peptide_template_anchors_xyz, anchor_tol, peptide_index, current_round)
	if(peptide_is_not_valid): return

	# 6. Create the peptide + MHC ensemble files
	peptide.create_peptide_receptor_complexes(filestore, receptor, peptide_index, current_round)

	# Done!
	return

def apegen(args):

	# 0. ARGUMENTS:
	
	parser = argparser.APE_Gen_parser()
	args = parser.parse_args()

	# - peptide_input: Crystal structure OR sequence
	peptide_input = args.peptide_input[0]

	# - receptor_class: .pdb OR sequence OR if peptide_input is crystal structure, REDOCK!
	receptor_class = args.receptor_class[0]

	# - Number of cores
	num_cores = int(args.num_cores)

	# - Number of loops on RCD (100 by default)
	num_loops = int(args.num_loops)

	# - RCD dist tolerance: RCD tolerance (in angstroms) of inner residues when performing IK
	RCD_dist_tol = args.RCD_dist_tol

	# - rigid_receptor : Disable sampling of receptor DoFs in the flex_res.txt
	doReceptorMinimization = not args.rigid_receptor

	# - Debug: Print extra information?
	debug = args.debug

	# --Save_only_pep_confs: Disable saving full conformations (peptide and MHC)
	saveFullConfs = not args.save_only_pep_confs

	# --Anchor_tolerance? Should this be an option?
	anchor_tol = args.anchor_tol

	# --Score_with_open_mm?
	score_with_openmm = args.score_with_openmm

	# - Number of rounds : When using multiple rounds, pass best scoring conformation across different rounds 
	num_rounds = args.num_rounds

	# (choose either 'receptor_only' or 'pep_and_recept')
	pass_type = args.pass_type

	# - min_with_smina: Minimize with SMINA instead of default Vinardo
	min_with_smina = args.min_with_smina

	# - use_gpu for Open_MM_minimization step
	device = "OpenCL" if args.use_gpu else "CPU"

	# --clean_rcd: Remove RCD folder at the end of each round?
	cleanRCD = args.clean_rcd

	# --force_restart: Force restart of APE-Gen *ONLY* in the first round and *ONLY* if no conformations are produced
	force_restart = args.force_restart

	# Directory to store intermediate files
	temp_files_storage = args.dir
	initialize_dir(temp_files_storage)

    # 1. INPUT PROCESSING

    # 1a. Peptide
	if debug: print("Processing Peptide Input")
	if peptide_input.endswith(".pdb"):
		# Fetch peptide sequence from .pdb and use that .pdb as a template
		peptide = Peptide.frompdb(peptide_input)
	else:
		# Fetch template from peptide template list
		peptide = Peptide.fromsequence(peptide_input)

	# Peptide Template is also a pMHC complex though	
	peptide_template = pMHC(pdb_filename = peptide.pdb_filename, peptide = peptide) 

	# The reason that we calculate the PTMs list here and not in the Peptide class is because we need it to be a
	# global variable to pass it on all peptide instances that need to be PTMed
	PTM_list = sequence_PTM_processing(peptide.sequence)
	peptide.sequence = re.sub('[a-z]', '', peptide.sequence)

	if debug: 
		print("Peptide Successfully Processed")
		print("Peptide Sequence: " + peptide.sequence)
		print("Peptide Template: " + peptide_template.pdb_filename)
		print("Peptide PTMs:")
		print(PTM_list)

	# 1b. Receptor
	if debug: print("Processing Receptor Input")
	if receptor_class.endswith(".pdb"):
		# If the file is .pdb, this will be your template! ##MUST CHECK VALIDITY IN THE FUNCTION
		receptor = Receptor.frompdb(receptor_class)
		receptor_template_file = receptor_class
	elif receptor_class.endswith(".fasta"):
		# If this is a sequence, the template is taken by MODELLER
		receptor = Receptor.fromfasta(receptor_class)
		receptor_template_file = receptor.pdb_filename
	elif receptor_class == "REDOCK":
		# If REDOCK, the receptor template is the peptide template!
		receptor = Receptor.fromredock(peptide_input)
		receptor_template_file = peptide.pdb_filename
	else:
		# If this is an allotype specification, fetch template like the peptide!
		receptor = Receptor.fromallotype(receptor_class)
		receptor_template_file = receptor.pdb_filename
	receptor.doMinimization = doReceptorMinimization
	receptor.useSMINA = min_with_smina

	# As receptor traditionally does not include the target, we will delete the .pdb for now:
	#receptor.pdb_filename = None (REMIND MYSELF WHY I DO THIS!)

	# Receptor Template is also a pMHC complex
	receptor_template = pMHC(pdb_filename = receptor_template_file, receptor = receptor) 
	if debug:
		print("Receptor Successfully Processed")
		print("Receptor Allotype: " + receptor.allotype)
		print("Receptor Template: " + receptor_template.pdb_filename)

	# 2. MAIN LOOP
	current_round = 1
	while current_round < num_rounds + 1: 

		# File storage location for the current round
		print("\n\nStarting round " + str(current_round) + " !!!!\n")
		filestore = temp_files_storage + "/" + str(current_round)

		# Alignment and preparing input for RCD
		# WARNING: pMHC complex at the time has both pMHC structures. 
		# It will be after the alignment that peptide file is a peptide and receptor file is a receptor
		if debug: print("Aligning peptide anchors to MHC pockets")
		receptor_template.align(reference = peptide_template, filestore = filestore)
		if debug: print("Preparing input to RCD")
		receptor_template.prepare_for_RCD(reference = peptide_template, filestore = filestore, pep_seq = peptide.sequence)
		receptor_template.add_sidechains(filestore = filestore)

		# Perform RCD on the receptor given peptide:
		if debug: print("Performing RCD")
		receptor_template.RCD(peptide, RCD_dist_tol, num_loops, filestore)


		# Prepare receptor for scoring (generate .pdbqt for SMINA):
		if debug: print("Preparing receptor for scoring (generate .pdbqt for SMINA)")
		initialize_dir(filestore + '/SMINA_data')
		receptor = receptor_template.receptor
		receptor.add_sidechains(filestore)
		receptor.prepare_for_scoring(filestore)

		# Get peptide template anchor positions for anchor tolerance filtering
		if debug: print("Extract peptide template anchors for anchor tolerance filtering")
		peptide_template_anchors_xyz = peptide_template.set_anchor_xyz(reference = peptide_template, 
																	   pep_seq = peptide.sequence)

		# Peptide refinement and scoring with SMINA on the receptor (done in parallel)
		if debug: print("Performing peptide refinement and scoring. This may take a while...")

		initialize_dir(filestore + '/SMINA_data/assembled_peptides') # Need to make initialize_dirs accepting a list
		initialize_dir(filestore + '/SMINA_data/per_peptide_results')
		initialize_dir(filestore + '/SMINA_data/PTMed_peptides')
		initialize_dir(filestore + '/SMINA_data/add_sidechains')
		initialize_dir(filestore + '/SMINA_data/pdbqt_peptides') 
		initialize_dir(filestore + '/SMINA_data/Scoring_results')
		initialize_dir(filestore + '/SMINA_data/flexible_receptors')
		initialize_dir(filestore + '/SMINA_data/minimized_receptors')
		initialize_dir(filestore + '/SMINA_data/Anchor_filtering')
		initialize_dir(filestore + '/SMINA_data/pMHC_complexes/')

		arg_list = list(map(lambda e: (e, filestore, PTM_list, receptor, peptide_template_anchors_xyz, anchor_tol, current_round), 
						range(1, num_loops + 1)))
		with WorkerPool(n_jobs=num_cores) as pool:
			results = pool.map(peptide_refinement_and_scoring, arg_list, progress_bar=True)

		initialize_dir(filestore + '/Final_conformations/')

		# Working Code for enbsembling individual .pdbs (might come in handy later!)

		# list_of_RCD_files = glob.glob(filestore + '/SMINA_data/PTMed_peptides/*.pdb')
		# ensembled = pdb_mkensemble.run(list_of_RCD_files)
		# with open(filestore + '/SMINA_data/PTMed_peptides/ensemble.pdb', 'w') as anchored_MHC_file:
		#	anchored_MHC_file.write(''.join(ensembled))

		# Code for non-parallel execution 

		#for argument in arg_list:
			#print(argument)
		#	peptide_refinement_and_scoring(argument[0], argument[1], argument[2], argument[3], argument[4], argument[5], argument[6])

		# Print and keep statistics
		print("\n\nEnd of round no. " + str(current_round) + "!!!")
		create_csv_from_list_of_files(filestore + '/total_results.csv', glob.glob(filestore + '/SMINA_data/per_peptide_results/*.log'))
		results_csv = pretty_print_analytics(filestore + '/total_results.csv')
		results_csv.to_csv(filestore + '/successful_conformations_statistics.csv', index = False)

		# OpenMM step. It seems to be parallel on its own, so no need to put it in a loop (but also inverstigate more?)		
		# A) Use nested for loops, as for every remaining conformation, we do *N* different OpenMM tries
		# B) Re-print results after the OpenMM step (include SMINA but also energies calculated from OpenMM?)
		if score_with_openmm:
			print("\n\n Chose OPENMM optimization!\n")
			initialize_dir(filestore + '/OpenMM_confs')
			initialize_dir(filestore + '/OpenMM_confs/minimized_receptors')
			initialize_dir(filestore + '/OpenMM_confs/minimized_complexes')
			initialize_dir(filestore + '/OpenMM_confs/pMHC_complexes/')

			successful_confs = results_csv['Peptide index'].tolist()
			for conf in tqdm(successful_confs, desc="pMHC conf", position = 0):
				
				# First I need to fix the receptors through PDBFixer:

				fixer = PDBFixer(filename=filestore + '/SMINA_data/minimized_receptors/receptor_' + str(conf) + ".pdb")
				fixer.findMissingResidues()
				fixer.removeHeterogens(True) #  True keeps water molecules while removing all other heterogens, REVISIT!
				fixer.findMissingAtoms()
				fixer.addMissingAtoms()
				fixer.addMissingHydrogens(7.0) # Ask Mauricio about those
				PDBFile.writeFile(fixer.topology, fixer.positions, open(filestore + '/OpenMM_confs/minimized_receptors/receptor_' + str(conf) + ".pdb", 'w'))

				# Unify peptide and receptor together
				pMHC_conformation = filestore + "/OpenMM_confs/pMHC_complexes/pMHC_" + str(conf) + ".pdb"
				merge_and_tidy_pdb([filestore + '/OpenMM_confs/minimized_receptors/receptor_' + str(conf) + ".pdb",
									filestore + '/SMINA_data/Anchor_filtering/peptide_' + str(conf) + ".pdb"], 
									pMHC_conformation)

				numTries = 10
				best_energy = float("inf")
				pMHC_complex = pMHC(pdb_filename = pMHC_conformation, peptide = peptide)
				for minimization_effort in tqdm(range(1, numTries + 1),  desc="No. of tries", position=1,
												leave=False):
					pMHC_complex = pMHC(pdb_filename = pMHC_conformation, peptide = peptide)
					best_energy = pMHC_complex.minimizeConf(filestore, best_energy, device)
			copy_batch_of_files(filestore + '/OpenMM_confs/minimized_complexes/', 
								filestore + '/Final_conformations/',
								query="pMHC_")
		else:
			copy_batch_of_files(filestore + '/SMINA_data/pMHC_complexes/', 
								filestore + '/Final_conformations',
								query="pMHC_")

		# Control whether there are no conformations. If they do, store the best one and continue.
		# If not, either abort or force restart (for round one)
		if(results_csv.shape[0] == 0 and current_round > 1):
			print('No conformations were produced this round. See results from previous rounds for final conformations. Aborting...')
			sys.exit(0)
		elif(results_csv.shape[0] == 0 and current_round == 1):
			if not force_restart:
				print('No conformations were produced... Aborting...')
				sys.exit(0)
			print('No conformations were produced... Force restarting...')
		else:
			# Storing the best conformation
			best_energy = results_csv['Affinity'].astype('float').min()
			best_conformation = results_csv[results_csv['Affinity'].astype('float') == best_energy]
			best_conformation_index = best_conformation['Peptide index'].values[0]
			print("\nStoring best conformation no. " + str(best_conformation_index) + " with Affinity = " + str(best_energy))
			copy_file(filestore + '/SMINA_data/pMHC_complexes/pMHC_' + str(best_conformation_index) + '.pdb',
				  	  filestore + '/min_energy_system.pdb')
			copy_file(filestore + '/SMINA_data/per_peptide_results/peptide_' + str(best_conformation_index) + '.log',
				  	  filestore + '/min_energy.log')

			# Decide where and how to pass the information for the next round
			receptor_template.pdb_filename = filestore + '/min_energy_system.pdb'
			if pass_type == 'pep_and_recept': peptide_template.pdb_filename = filestore + '/min_energy_system.pdb'

			# Finally, advance to the next round!
			current_round += 1
			
	# Ending and final statisticsI can 
	print("\n\nEnd of APE-Gen !!!")
	create_csv_from_list_of_files(temp_files_storage + '/APE_gen_best_run_results.csv', 
								  ["{}/{}/min_energy.log".format(temp_files_storage,i) for i in range(1, num_rounds + 1)])
	results_csv = pretty_print_analytics(temp_files_storage + '/APE_gen_best_run_results.csv')

if __name__ == "__main__":
    apegen(sys.argv[1:])