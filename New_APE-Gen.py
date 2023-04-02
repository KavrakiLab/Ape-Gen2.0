from helper_scripts import argparser
from helper_scripts.Ape_gen_macros import apply_function_to_file, replace_chains, initialize_dir,  \
											merge_and_tidy_pdb, add_sidechains,					   \
											create_csv_from_list_of_files, 						   \
											copy_file, pretty_print_analytics, move_batch_of_files,\
											copy_batch_of_files, split_receptor_and_peptide,	   \
											remove_remarks_and_others_from_pdb, replace_HETATM,    \
											delete_elements, verbose, set_verbose                  \

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

def rescoring_after_openmm(conf_index, filestore, rcd_num_loops, current_round, peptide_template_anchors_xyz, anchor_tol, tolerance_anchors):

	new_filestore = filestore + '/5_openMM_conformations'

	# 1. Rename B chain to C chain
	apply_function_to_file(replace_chains, new_filestore + "/10_pMHC_complexes/pMHC_" + str(conf_index) + ".pdb", chain_from="B", chain_to="C")

	# 2. Separate the peptide from the MHC
	receptor_file, peptide_file = split_receptor_and_peptide(new_filestore + "/10_pMHC_complexes/pMHC_" + str(conf_index) + ".pdb")

	# 3. Prepare Receptor for Scoring
	apply_function_to_file(remove_remarks_and_others_from_pdb, receptor_file, records=('ATOM', 'HETATM', 'TER', 'END '))

	receptor = Receptor.frompdb(receptor_file)
	receptor.doMinimization = True
	receptor.useSMINA = True
	receptor.prepare_for_scoring(new_filestore + '/09_minimized_receptors', index=str(conf_index))

	apply_function_to_file(remove_remarks_and_others_from_pdb, peptide_file, records=('ATOM', 'HETATM', 'TER', 'END '))

	apply_function_to_file(replace_HETATM, peptide_file)
	
	peptide = Peptide.frompdb(peptide_file, secondary_anchors = tolerance_anchors, peptide_index=conf_index)
	peptide.sequence = re.sub('[a-z]', '', peptide.sequence) # Remove PTMs from the sequence

	# 4. Re-score with SMINA (enforce no further minimization)
	peptide_is_not_valid = peptide.prepare_for_scoring(new_filestore, current_round)
	if(peptide_is_not_valid): return

	peptide.score_with_SMINA(new_filestore, receptor)

	# 5. Anchor filtering step (probably not needed, anchors are not moving that much)
	peptide_is_not_valid = peptide.compute_anchor_tolerance(new_filestore, receptor, peptide_template_anchors_xyz, anchor_tol, current_round, rcd_num_loops)
	if(peptide_is_not_valid): return

	# 6. Create the peptide + MHC ensemble files (Already have those but ok...)
	peptide.create_peptide_receptor_complexes(new_filestore, receptor)

	# Done!
	return

def prepare_for_openmm(conf_index, filestore, peptide):

	# 1. Run receptors through PDBFixer, as the non-polar hydrogens could be in wrong places:
	receptor = Receptor.frompdb(filestore + '/4_SMINA_data/09_minimized_receptors/receptor_' + str(conf_index) + ".pdb")
	add_sidechains(receptor.pdb_filename, filestore, add_hydrogens="Yes", keep_IDs=True)

	# 2. Unify peptide and receptor together and create a new pMHC complex
	pMHC_conformation = filestore + "/5_openMM_conformations/11_pMHC_before_sim/pMHC_" + str(conf_index) + ".pdb"
	merge_and_tidy_pdb([receptor.pdb_filename,
						filestore + '/4_SMINA_data/08_anchor_filtering/peptide_' + str(conf_index) + ".pdb"],
						pMHC_conformation)
	pMHC_complex = pMHC(pdb_filename=pMHC_conformation, peptide=peptide)

	# 3. If there is a phosphorylation somewhere, we need to give the appropriate CONECT fields to the PTM residue
	pMHC_complex.add_PTM_CONECT_fields(filestore, peptide.PTM_list, conf_index)

	# Done!
	return

def peptide_refinement_and_scoring(index, rcd_num_loops, original_peptide, filestore, receptor, tolerance_anchors, peptide_template_anchors_xyz, anchor_tol, current_round):

	# Routine that refines and scores a peptide/receptor pair with SMINA/Vinardo
	new_filestore = filestore + '/4_SMINA_data'

	# 1. Assemble peptide by mergin the peptide anchors and the middle part
	assembled_peptide = new_filestore + '/01_assembled_peptides/assembled_' + str(index) + '.pdb'
	if index <= rcd_num_loops:
		model_location = filestore + '/3_RCD_data/splits/model_' + str(index) + '.pdb'
		Nterm_location = filestore + '/2_input_to_RCD/N_ter.pdb'
		Cterm_location = filestore + '/2_input_to_RCD/C_ter.pdb'
		merge_and_tidy_pdb([Nterm_location, model_location, Cterm_location], assembled_peptide)
	else:
		copy_file(filestore + '/2_input_to_RCD/model_0.pdb', assembled_peptide)

	peptide = Peptide.frompdb(assembled_peptide, secondary_anchors=tolerance_anchors, peptide_index=index, PTM_list = original_peptide.PTM_list)

	# 2. Now that the peptide is assembled, Fill in the sidechains with pdbfixer
	peptide.pdb_filename = add_sidechains(peptide.pdb_filename, new_filestore, peptide_idx=index)

	# 3. Do PTMs
	peptide_is_not_valid = peptide.perform_PTM(new_filestore, current_round)
	if(peptide_is_not_valid): return

	# 4. Score with SMINA
	# 4a. .pdb to .pdbqt transformation using autodocktools routines (very good for filtering bad conformations)
	peptide_is_not_valid = peptide.prepare_for_scoring(new_filestore, current_round)
	if(peptide_is_not_valid): return

	# 4b. Optimize and score with SMINA (or other options, depending on args)
	peptide.dock_score_with_SMINA(new_filestore, receptor)

	# 5. Anchor filtering (based on anchor tolerance argument) (improvement from previous version)
	peptide_is_not_valid = peptide.compute_anchor_tolerance(new_filestore, receptor, peptide_template_anchors_xyz, anchor_tol, current_round, rcd_num_loops)
	if(peptide_is_not_valid): return

	# 6. Fix flexible residue co-ordinates if receptor is flexible
	if receptor.doMinimization:
		peptide_is_not_valid = peptide.fix_flexible_residues(new_filestore, receptor, current_round)
		if(peptide_is_not_valid): return

	# 7. Create the peptide + MHC ensemble files
	peptide.create_peptide_receptor_complexes(new_filestore, receptor)

	# Done!
	return

def apegen(args):
	print("Start of APE-Gen")

	# 0. ARGUMENTS:

	parser = argparser.APE_Gen_parser()
	args = parser.parse_args()

	# - peptide_input: Crystal structure OR sequence
	peptide_input = args.peptide_input[0]

	# - receptor_class: .pdb OR sequence OR if peptide_input is crystal structure, REDOCK!
	receptor_class = args.receptor_class[0]

	# - Number of cores
	num_cores = int(args.num_cores)

	# - Number of loops to generate with RCD
	rcd_num_loops = int(args.num_generated_loops)

	# - Number of loops to optimize (that will pass as a result of a loop scoring function)
	num_loops = int(args.num_loops_for_optimization)

	# The percentage of overall peptide conformations processed (defined by --num_loops_for_optimization flag) that will be coming from RCD sampling.
	non_sampled_confs = int(np.rint((1 - args.sampling_ratio)*num_loops))
	num_loops = int(np.rint(args.sampling_ratio*num_loops))

	# - RCD dist tolerance: RCD tolerance (in angstroms) of inner residues when performing IK
	RCD_dist_tol = args.RCD_dist_tol

	# --loop_score: Choose scoring function for RCD loop scoring (none will avoid scoring altogether)
	loop_score = args.loop_score

	# - rigid_receptor : Disable sampling of receptor DoFs in the flex_res.txt
	doReceptorMinimization = not args.rigid_receptor

	# - Debug: Print extra information?
	verbose = args.verbose
	set_verbose(verbose)

	# --Save_only_pep_confs: Disable saving full conformations (peptide and MHC)
	saveFullConfs = not args.save_only_pep_confs

	# --anchors: User defined anchors for peptide template search + anchor tolerance
	anchors = args.anchors

	# --Anchor_tolerance? Should this be an option?
	anchor_tol = args.anchor_tol

	# --Score_with_open_mm?
	score_with_openmm = args.score_with_openmm

	# Do not apply constraints on the backbone when applying openMM
	no_constraints_openmm = args.no_constraints_openmm

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

	# --anchor_selection: Give what type of anchors should be considered in the anchor tolerance step (choose 'primary', 'secondary' or 'none' to skip the anchor tolerance step altogether)
	anchor_selection = args.anchor_selection

	# --cv: ONLY FOR TESTING (to be removed in the final version)
	cv = args.cv

	# Directory to store intermediate files
	temp_files_storage = args.dir
	initialize_dir(temp_files_storage)

	# 1. INPUT PROCESSING
	peptide = Peptide.init_peptide(peptide_input)
	PTM_list = peptide.PTM_list

	receptor, receptor_template_file = Receptor.init_receptor(receptor_class, temp_files_storage +  '/MODELLER_output', peptide.sequence, cv)
	receptor.doMinimization = doReceptorMinimization
	receptor.useSMINA = min_with_smina
	
	# Peptide Template and Receptor Template are pMHC complexes	
	peptide_template = peptide.get_peptide_template(receptor.allotype, anchors, anchor_selection, cv)
	receptor_template = pMHC(pdb_filename=receptor_template_file, peptide=peptide, receptor=receptor)

	# Get peptide template anchor positions for anchor tolerance filtering
	peptide_template_anchors_xyz, tolerance_anchors = peptide_template.set_anchor_xyz(anchor_selection, peptide)

	if verbose:
		print("\nReceptor Successfully Processed")
		print("    Receptor Allotype: " + receptor.allotype)
		print("    Receptor Template: " + receptor_template.pdb_filename)

		print("\nPeptide Successfully Processed")
		print("    Peptide Sequence: " + peptide.sequence)
		print("    Peptide Template: " + peptide_template.pdb_filename)
		print("    Peptide Primary Anchors:")
		print("       ", peptide.primary_anchors)
		print("    Peptide Secondary Anchors:")
		print("       ", peptide.secondary_anchors)
		print("    Peptide PTMs:")
		print("       ", PTM_list)

	# Check if:
	# A. There are any PTMs other than phosphorylation. GROMACS will be considered, but not right know..
	# B. Phosphorylation is on N-terminus or C-terminus. FF parameters are not given for these cases.
	# C. User wants to model with no hydrogens involved, but also run an energy minimization routine. 
	#    From my understanding, PDBFixer when given an MHC with no hydrogens will mess smth up not in terms of atoms, but in terms of bonds. 
	#    Let's prevent users from actually doing this.
	# D. If the number of loops is more than 5000, we should remove possibility from modelling (and this I guess can be manually overriden) 
	if (('phosphorylate 1' in PTM_list) or ('phosphorylate ' + str(len(peptide.sequence)) in PTM_list)) and (score_with_openmm):
		sys.exit("\nERROR: Phosphorylation in N-terminus or C-terminus and openMM optimization is NOT supported. Force Field parameters are not released yet. Please omit OpenMM step for modelling this type of PTM.")
	for PTM in PTM_list:
		if (not PTM.startswith("phosphorylate")) and (score_with_openmm):
			sys.exit("\nERROR: PTM other than phosphorylation is not yet supported with OpenMM. Omit the OpenMM step and stay tuned for changes!")
	if num_loops > rcd_num_loops:
		sys.exit("\nERROR: The number of loops for post-processing should not exceed the number of loops that are generated!")

	# 2. MAIN LOOP
	current_round = 1
	while current_round < num_rounds + 1:

		# File storage location for the current round
		if verbose: print("\n\nStarting round " + str(current_round) + " !!!!\n")
		filestore = temp_files_storage + "/" + str(current_round)

		# Alignment and preparing input for RCD
		# WARNING: pMHC complex at the time has both pMHC structures.
		# It will be after the alignment that peptide file is a peptide and receptor file is a receptor
		if verbose: print("Aligning peptide anchors to MHC pockets")
		receptor_template.align(reference=peptide_template, filestore=filestore)
		if verbose: print("Preparing input to RCD")
		receptor_template.prepare_for_RCD_v3(reference=peptide_template, peptide=peptide, filestore=filestore)
		add_sidechains(receptor_template.pdb_filename, filestore, keep_IDs=True)

		# Perform RCD on the receptor given peptide:
		if verbose: print("Performing RCD")
		loop_indexes = receptor_template.RCD(peptide, RCD_dist_tol, rcd_num_loops, num_loops, 
											 loop_score, non_sampled_confs, filestore)

		# Prepare receptor for scoring (generate .pdbqt for SMINA):
		if verbose: print("Preparing receptor for scoring (generate .pdbqt for SMINA)")
		initialize_dir(filestore + '/4_SMINA_data')
		receptor = receptor_template.receptor
		add_sidechains(receptor.pdb_filename, filestore)
		receptor.prepare_for_scoring(filestore + "/4_SMINA_data")

		# Peptide refinement and scoring with SMINA on the receptor (done in parallel)
		if verbose: print("Performing peptide refinement and scoring. This may take a while...")

		subdir_list = ['/01_assembled_peptides', '/05_per_peptide_results', '/03_PTMed_peptides',
						'/02_add_sidechains', '/04_pdbqt_peptides', '/06_scoring_results', '/07_flexible_receptors',
						'/09_minimized_receptors', '/08_anchor_filtering', '/10_pMHC_complexes/']
		initialize_dir([filestore + '/4_SMINA_data' + subdir for subdir in subdir_list])
		
		arg_list = list(map(lambda pep_index: (pep_index, rcd_num_loops, peptide, filestore, receptor, tolerance_anchors, peptide_template_anchors_xyz, anchor_tol, current_round), 
						loop_indexes))
		with WorkerPool(n_jobs=num_cores) as pool:
			results = pool.map(peptide_refinement_and_scoring, arg_list, progress_bar=verbose)

		# Code for non-parallel execution and debugging

		#for argument in arg_list:
		#    print(argument)
		#    peptide_refinement_and_scoring(argument[0], argument[1], argument[2], argument[3], argument[4], argument[5], argument[6], argument[7], argument[8])

		# Print and keep statistics
		best_conf_dir = filestore + '/4_SMINA_data'
		if verbose: print("\n\nEnd of main workflow of round no. " + str(current_round) + "!!!")
		create_csv_from_list_of_files(filestore + '/4_SMINA_data/total_results.csv', glob.glob(filestore + '/4_SMINA_data/05_per_peptide_results/*.log'))
		results_csv = pretty_print_analytics(filestore + '/4_SMINA_data/total_results.csv', verbose=verbose)
		results_csv.to_csv(filestore + '/4_SMINA_data/successful_conformations_statistics.csv', index=False)

		# OpenMM step
		if(score_with_openmm and results_csv.shape[0] > 0):

			if verbose: print("\n\nOpennMM optimization!\n")

			dir_list = ['/6_final_conformations/', '/5_openMM_conformations']
			subdir_list = ['/fixed_receptors', '/14_minimized_complexes', '/10_pMHC_complexes', 
							'/11_pMHC_before_sim', '/13_connected_pMHC_complexes', '/12_PTM_conect_indexes',
							'/06_scoring_results', '/05_per_peptide_results', '/09_minimized_receptors', 
							'/08_anchor_filtering', '/04_pdbqt_peptides', '/07_flexible_receptors']

			initialize_dir([filestore + dir for dir in dir_list])
			initialize_dir([filestore + '/5_openMM_conformations' + subdir for subdir in subdir_list])

			successful_confs = results_csv['Peptide index'].tolist()
			if verbose: print("Preparing input for OpenMM optimization. This may take a while...")

			# First prepare for OpenMM
			arg_list = list(map(lambda e: (e, filestore, peptide), successful_confs))
			with WorkerPool(n_jobs=min(num_cores, len(successful_confs))) as pool:
				results = pool.map(prepare_for_openmm, arg_list, progress_bar=verbose)

			# Actual minimization step
			if verbose:
				print("\nMinimizing energy...")
				if no_constraints_openmm: print("Removing backbone constraints from energy minimization!")
				disable_progress_bar = False
				leave_progress_bar=False
			else:
				disable_progress_bar = True
				leave_progress_bar=True

			for conf_index in tqdm(successful_confs, desc="pMHC conf", position=0, disable=disable_progress_bar):

				numTries = 1
				best_energy = float("inf")
				pMHC_complex = pMHC(pdb_filename=filestore + "/5_openMM_conformations/13_connected_pMHC_complexes/pMHC_" + str(conf_index) + ".pdb", 
									peptide=peptide)
				for minimization_effort in tqdm(range(1, numTries + 1),  desc="No. of tries", position=1,
												leave=leave_progress_bar, disable=disable_progress_bar):
					best_energy = pMHC_complex.minimizeConf(filestore, best_energy, no_constraints_openmm, device)
					with open(filestore + "/5_openMM_conformations/05_per_peptide_results/peptide_" + str(conf_index) + ".log", 'w') as peptide_handler:
						peptide_handler.write(str(current_round) + "," + str(conf_index) + ",Successfully Modeled," + str(best_energy) + "\n")

			# Rescoring and re-filtering resulting conformations
			if verbose: print("\nRescoring and re-filtering resulting conformations:")
			arg_list = list(map(lambda conf_index: (conf_index, filestore, rcd_num_loops, current_round, peptide_template_anchors_xyz, anchor_tol, tolerance_anchors), successful_confs))
			with WorkerPool(n_jobs=min(num_cores, len(successful_confs))) as pool:
				results = pool.map(rescoring_after_openmm, arg_list, progress_bar=verbose)

			copy_batch_of_files(filestore + '/5_openMM_conformations/10_pMHC_complexes/',
								filestore + '/6_final_conformations/',
								query="pMHC_")

			best_conf_dir = filestore + '/5_openMM_conformations'
			if verbose: print("\n\nEnd of OpenMM step of round no. " + str(current_round) + "!!!")
			create_csv_from_list_of_files(filestore + '/5_openMM_conformations/total_results.csv', glob.glob(filestore + '/5_openMM_conformations/05_per_peptide_results/*.log'))
			results_csv = pretty_print_analytics(filestore + '/5_openMM_conformations/total_results.csv', verbose=verbose)
			results_csv.to_csv(filestore + '/5_openMM_conformations/successful_conformations_statistics.csv', index=False)

		else:
			initialize_dir(filestore + '/5_final_conformations/')
			copy_batch_of_files(filestore + '/4_SMINA_data/10_pMHC_complexes/',
								filestore + '/5_final_conformations',
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
			anchor_tol += 1
			if anchor_tol > 10:
				print('No conformations were produced... Tolerance is quite high so something else is wrong... Aborting...')
				sys.exit(0)
			print('No conformations were produced... Force restarting with increased anchor tolerance = ' + str(anchor_tol))
			receptor_template.pdb_filename = receptor_template_file
		else:
			# Storing the best conformation
			best_energy = results_csv['Affinity'].astype('float').min()
			best_conformation = results_csv[results_csv['Affinity'].astype('float') == best_energy]
			best_conformation_index = best_conformation['Peptide index'].values[0]
			if verbose: print("\nStoring best conformation no. " + str(best_conformation_index) + " with Affinity = " + str(best_energy))
			copy_file(best_conf_dir + '/10_pMHC_complexes/pMHC_' + str(best_conformation_index) + '.pdb',
					  best_conf_dir + '/min_energy_system.pdb')
			copy_file(best_conf_dir + '/05_per_peptide_results/peptide_' + str(best_conformation_index) + '.log',
					  best_conf_dir + '/min_energy.log')

			# Decide where and how to pass the information for the next round
			receptor_template.pdb_filename = best_conf_dir + '/min_energy_system.pdb'
			if pass_type == 'pep_and_recept': peptide_template.pdb_filename = best_conf_dir + '/min_energy_system.pdb'

			# Finally, advance to the next round!
			current_round += 1

	# Ending and final statistics
	print("\n\nEnd of APE-Gen")
	best_conf_dir = '5_openMM_conformations' if score_with_openmm else '4_SMINA_data'
	create_csv_from_list_of_files(temp_files_storage + '/APE_gen_best_run_results.csv',
								  ["{}/{}/{}/min_energy.log".format(temp_files_storage,i,best_conf_dir) for i in range(1, num_rounds + 1)])
	results_csv = pretty_print_analytics(temp_files_storage + '/APE_gen_best_run_results.csv', verbose=True)

if __name__ == "__main__":
	apegen(sys.argv[1:])
