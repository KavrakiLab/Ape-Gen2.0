from helper_scripts import argparser
from helper_scripts.Ape_gen_macros import apply_function_to_file, replace_chains, initialize_dir,  \
											merge_and_tidy_pdb, add_sidechains,					   \
											create_csv_from_list_of_files, 						   \
											copy_file, pretty_print_analytics, move_batch_of_files,\
											copy_batch_of_files, split_receptor_and_peptide,	   \
											remove_remarks_and_others_from_pdb, replace_HETATM,    \
											delete_elements, split_to_equal_parts, remove_dirs,    \
											verbose, set_verbose                  \

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

def rescoring_after_openmm(conf_index, filestore, rcd_num_loops, peptide_template_anchors_xyz, anchor_tol, tolerance_anchors, min_with_smina):

	new_filestore = filestore + '/5_openMM_conformations'

	# 1. Rename B chain to C chain
	apply_function_to_file(replace_chains, new_filestore + "/10_pMHC_complexes/pMHC_" + conf_index + ".pdb", chain_from="B", chain_to="C")

	# 2. Separate the peptide from the MHC
	receptor_file, peptide_file = split_receptor_and_peptide(new_filestore + "/10_pMHC_complexes/pMHC_" + conf_index + ".pdb")

	# 3. Prepare Receptor for Scoring
	apply_function_to_file(remove_remarks_and_others_from_pdb, receptor_file, records=('ATOM', 'HETATM', 'TER', 'END '))

	receptor = Receptor.frompdb(receptor_file)
	receptor.doMinimization = True
	receptor.useSMINA = min_with_smina
	receptor_is_not_valid = receptor.prepare_for_scoring(new_filestore + '/09_minimized_receptors', index=conf_index)
	if(receptor_is_not_valid): return

	apply_function_to_file(remove_remarks_and_others_from_pdb, peptide_file, records=('ATOM', 'HETATM', 'TER', 'END '))

	apply_function_to_file(replace_HETATM, peptide_file)
	
	peptide = Peptide.frompdb(peptide_file, secondary_anchors = tolerance_anchors, peptide_index=conf_index)
	peptide.sequence = re.sub('[a-z]', '', peptide.sequence) # Remove PTMs from the sequence

	# 4. Re-score with SMINA (enforce no further minimization)
	peptide_is_not_valid = peptide.prepare_for_scoring(new_filestore)
	if(peptide_is_not_valid): return

	peptide.score_with_SMINA(new_filestore, receptor)

	# 5. Anchor filtering step (probably not needed, anchors are not moving that much)
	peptide_is_not_valid = peptide.compute_anchor_tolerance(new_filestore, receptor, peptide_template_anchors_xyz, anchor_tol, rcd_num_loops)
	if(peptide_is_not_valid): return

	# 6. Create the peptide + MHC ensemble files (Already have those but ok...)
	peptide.create_peptide_receptor_complexes(new_filestore, receptor)

	# Done!
	return

def prepare_for_openmm(conf_index, filestore, peptide):

	# 1. Run receptors through PDBFixer, as the non-polar hydrogens could be in wrong places:
	receptor = Receptor.frompdb(filestore + '/4_SMINA_data/09_minimized_receptors/receptor_' + conf_index + ".pdb")
	add_sidechains(receptor.pdb_filename, filestore, add_hydrogens="Yes", keep_IDs=True)

	# 2. Unify peptide and receptor together and create a new pMHC complex
	pMHC_conformation = filestore + "/5_openMM_conformations/11_pMHC_before_sim/pMHC_" + conf_index + ".pdb"
	merge_and_tidy_pdb([receptor.pdb_filename,
						filestore + '/4_SMINA_data/08_anchor_filtering/peptide_' + conf_index + ".pdb"],
						pMHC_conformation)
	pMHC_complex = pMHC(pdb_filename=pMHC_conformation, peptide=peptide)

	# 3. If there is a phosphorylation somewhere, we need to give the appropriate CONECT fields to the PTM residue
	pMHC_complex.add_PTM_CONECT_fields(filestore, peptide.PTM_list, conf_index)

	# Done!
	return

def peptide_refinement_and_scoring(index, template_index, new_index, rcd_num_loops, original_peptide, filestore, receptor, tolerance_anchors, peptide_template_anchors_xyz, anchor_tol):

	# Routine that refines and scores a peptide/receptor pair with SMINA/Vinardo
	new_filestore = filestore + '/4_SMINA_data'

	# 1. Assemble peptide by mergin the peptide anchors and the middle part
	assembled_peptide = new_filestore + '/01_assembled_peptides/assembled_' + new_index + '.pdb'
	if index <= rcd_num_loops:
		model_location = filestore + '/3_RCD_data/' + str(template_index) + '/splits/model_' + str(index) + '.pdb'
		Nterm_location = filestore + '/2_input_to_RCD/' + str(template_index) + '/N_ter.pdb'
		Cterm_location = filestore + '/2_input_to_RCD/' + str(template_index) + '/C_ter.pdb'
		merge_and_tidy_pdb([Nterm_location, model_location, Cterm_location], assembled_peptide)
	else:
		copy_file(filestore + '/2_input_to_RCD/' + str(template_index) + '/model_'+ str(template_index) + '.pdb',
				  assembled_peptide)

	peptide = Peptide.frompdb(assembled_peptide, secondary_anchors=tolerance_anchors, peptide_index=new_index, PTM_list = original_peptide.PTM_list)

	# 2. Now that the peptide is assembled, Fill in the sidechains with pdbfixer
	peptide.pdb_filename = add_sidechains(peptide.pdb_filename, new_filestore, peptide_idx=new_index, keep_IDs=True)

	# 3. Do PTMs
	peptide_is_not_valid = peptide.perform_PTM(new_filestore)
	if(peptide_is_not_valid): return

	# 4. Score with SMINA
	# 4a. .pdb to .pdbqt transformation using autodocktools routines (very good for filtering bad conformations)
	peptide_is_not_valid = peptide.prepare_for_scoring(new_filestore)
	if(peptide_is_not_valid): return

	# 4b. Optimize and score with SMINA (or other options, depending on args)
	peptide_is_not_valid = peptide.dock_score_with_SMINA(new_filestore, receptor)
	if(peptide_is_not_valid): return

	# 5. Anchor filtering (based on anchor tolerance argument) (improvement from previous version)
	peptide_is_not_valid = peptide.compute_anchor_tolerance(new_filestore, receptor, peptide_template_anchors_xyz, anchor_tol, rcd_num_loops)
	if(peptide_is_not_valid): return

	# 6. Fix flexible residue co-ordinates if receptor is flexible
	if receptor.doMinimization:
		peptide_is_not_valid = peptide.fix_flexible_residues(new_filestore, receptor)
		if(peptide_is_not_valid): return

	# 7. Create the peptide + MHC ensemble files
	peptide.create_peptide_receptor_complexes(new_filestore, receptor)

	# Done!
	return

def backbone_sampling(template_index, peptide_templates, receptor_template, peptide, anchors, anchor_status, anchor_selection, rcd_num_loops, RCD_dist_tol, filestore):		
		
	# Routine that initializes a peptide template and samples backbones based no that template

	# 1. Initialize appropriate directories
	initialize_dir([filestore + '/1_alignment_files/' + str(template_index),
					filestore + '/2_input_to_RCD/' + str(template_index),
					filestore + '/3_RCD_data/'+ str(template_index)])
		
	# 2. Peptide template initialization
	peptide_template = peptide.initialize_peptide_template(peptide_templates.iloc[[template_index]], anchors, anchor_status)

	# 3. Alignment and preparing input for RCD
	receptor_template.align(reference=peptide_template, filestore=filestore, template_index=template_index)

	# 4. Get peptide template anchor positions for anchor tolerance filtering
	#anchor_filtering_data_dict[template_index] = peptide_template.set_anchor_xyz(anchor_selection, peptide)
	receptor_template.prepare_for_RCD(reference=peptide_template, peptide=peptide, 
											 filestore=filestore, template_index=template_index)

	# 5. Perform RCD on the receptor given peptide:
	receptor_template.RCD(peptide, RCD_dist_tol, rcd_num_loops, filestore, template_index)

	# Done! (6. Return also the anchor information for each template, it will come in handy later on)
	return (template_index, peptide_template.set_anchor_xyz(anchor_selection, peptide))

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

	# - rigid_receptor : Disable sampling of receptor DoFs in the ./helper_files/flex_res.txt
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

	# (choose either 'receptor_only' or 'pep_and_recept')
	pass_type = args.pass_type

	# - min_with_smina: Minimize with SMINA instead of default Vinardo
	min_with_smina = args.min_with_smina

	# - use_gpu for Open_MM_minimization step
	device = "OpenCL" if args.use_gpu else "CPU"

	# --clean_rcd: Remove RCD folder at the end of each round?
	cleanRCD = args.clean_rcd

	# --anchor_selection: Give what type of anchors should be considered in the anchor tolerance step (choose 'primary', 'secondary' or 'none' to skip the anchor tolerance step altogether)
	anchor_selection = args.anchor_selection

	# --max_no_templates: The maximum number of templates that will be used in the modelling process.
	max_no_templates = args.max_no_templates

	# --similarity_threshold: Score [0-1] that defines if a peptide template will be considered as a candidate during the modelling process.
	similarity_threshold = args.similarity_threshold

	# Option for anchor identification: Either PMBEC or MHCflurry motifs
	use_motifs = args.use_motifs

	# --keep_all_files: Keep all intermediate generated files from the modeling process
	keep_all_files = args.keep_all_files

	# --cv: ONLY FOR TESTING (to be removed in the final version)
	cv = args.cv

	# Directory to store intermediate files
	temp_files_storage = args.dir
	initialize_dir(temp_files_storage)

	# 1. INPUT PROCESSING
	peptide = Peptide.init_peptide(peptide_input)
	PTM_list = peptide.PTM_list

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

	# File storage location
	filestore = temp_files_storage + "/results"

	if verbose: print("Initializing receptor")
	receptor, receptor_template_file = Receptor.init_receptor(receptor_class, temp_files_storage +  '/MODELLER_output', peptide.sequence, cv)
	receptor.doMinimization = doReceptorMinimization
	receptor.useSMINA = min_with_smina
	
	# Peptide Template and Receptor Template are pMHC complexes	
	peptide_templates, anchors, anchor_status = peptide.get_peptide_templates(receptor.allotype, anchors, max_no_templates, similarity_threshold, use_motifs, cv)
	if peptide_templates.empty:
		print("No available peptides for the given peptide-MHC input! Aborting...")
		sys.exit(0)

	# Prepare receptor for scoring (generate .pdbqt for SMINA):
	receptor_template = pMHC(pdb_filename=receptor_template_file, peptide=peptide, receptor=receptor)

	if verbose:
		print("\nReceptor Successfully Processed")
		print("    Receptor Allotype: " + receptor.allotype)
		print("    Receptor Template: " + receptor_template.pdb_filename)
		print("\nPeptide Successfully Processed")
		print("    Peptide Sequence:")
		print("       ", peptide.sequence)
		print("    Peptide PTMs:")
		print("       ", PTM_list)

	
	anchor_filtering_data_dict = {}
	filestore = temp_files_storage + "/results/"
	
	# 2. BACKBONE SAMPLING LOOP
	arg_list = list(map(lambda template_index: (template_index, peptide_templates, receptor_template, peptide, anchors, anchor_status, anchor_selection, rcd_num_loops, RCD_dist_tol, filestore), 
						list(range(peptide_templates.shape[0]))))
	with WorkerPool(n_jobs=num_cores) as pool:
		anchor_results = pool.map(backbone_sampling, arg_list, progress_bar=verbose)

	# 3. LOOP SCORING LOOP
	if verbose: print("Scoring the sampled loops...")
	loop_index_list = []
	new_index_list = []
	template_index_list = []
	num_loops_list = split_to_equal_parts(num_loops, peptide_templates.shape[0])
	non_sampled_confs_list = split_to_equal_parts(non_sampled_confs, peptide_templates.shape[0])
	anchor_filtering_data_dict = {}

	for template_index in range(peptide_templates.shape[0]):
		
		# Extract anchor information from mpire process
		anchor_info = anchor_results[template_index]
		anchor_filtering_data_dict[anchor_info[0]] = (anchor_info[1][0], anchor_info[1][1])

		# Rank the loops and make indexes intepretable (template index + loop index)
		loop_indexes = receptor_template.loop_ranking(rcd_num_loops, num_loops_list[template_index], loop_score, non_sampled_confs_list[template_index], filestore, template_index)
		no_of_conformations = num_loops_list[template_index] + non_sampled_confs_list[template_index]
		new_indexes = [str(template_index) + str(i).zfill(len(str(no_of_conformations))) for i in range(non_sampled_confs_list[template_index], no_of_conformations)]
		new_indexes = new_indexes + [str(template_index) + str(i).zfill(len(str(no_of_conformations))) for i in range(0, non_sampled_confs_list[template_index])]
		loop_index_list += loop_indexes
		new_index_list += new_indexes
		template_index_list += [template_index]*no_of_conformations

	# 4. PEPTIDE REFINEMENT AND SCORING LOOP
	subdir_list = ['/01_assembled_peptides', '/05_per_peptide_results', '/03_PTMed_peptides',
				   '/02_add_sidechains', '/04_pdbqt_peptides', '/06_scoring_results', '/07_flexible_receptors',
				   '/09_minimized_receptors', '/08_anchor_filtering', '/10_pMHC_complexes/']
	initialize_dir([filestore + '/4_SMINA_data' + subdir for subdir in subdir_list])

	if verbose: print("Preparing receptor sans peptide for scoring (generate receptor.pdbqt)")
	receptor_template.remove_peptide(filestore + "/4_SMINA_data")
	receptor = receptor_template.receptor
	add_sidechains(receptor.pdb_filename, filestore, keep_IDs=True)
	receptor_is_not_valid = receptor.prepare_for_scoring(filestore + "/4_SMINA_data")
	if(receptor_is_not_valid):
		print("There is something wrong with the receptor file... Check the logs! Aborting...")
		sys.exit(0)
	
	if verbose: print("Performing peptide refinement and scoring. This may take a while...")
	arg_list = []
	for i, pep_index in enumerate(loop_index_list):
		arg_list.append((pep_index, template_index_list[i], new_index_list[i], rcd_num_loops, peptide, filestore, 
						receptor, anchor_filtering_data_dict[template_index_list[i]][1], 
						anchor_filtering_data_dict[template_index_list[i]][0], anchor_tol))
	with WorkerPool(n_jobs=num_cores) as pool:
		results = pool.map(peptide_refinement_and_scoring, arg_list, progress_bar=verbose)

	# Code for non-parallel execution and debugging
	#for argument in arg_list:
	#    print(argument)
	#    peptide_refinement_and_scoring(argument[0], argument[1], argument[2], argument[3], argument[4], argument[5], argument[6], argument[7], argument[8])

	# Print and keep statistics
	best_conf_dir = filestore + '/4_SMINA_data'
	if verbose: print("\n\nEnd of main workflow !!!")
	create_csv_from_list_of_files(filestore + '/4_SMINA_data/total_results.csv', glob.glob(filestore + '/4_SMINA_data/05_per_peptide_results/*.log'))
	results_csv = pretty_print_analytics(filestore + '/4_SMINA_data/total_results.csv', verbose=verbose)
	results_csv.to_csv(temp_files_storage + '/successful_conformations_statistics.csv', index=False)

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
			leave_progress_bar = False
		else:
			disable_progress_bar = True
			leave_progress_bar = True

		for conf_index in tqdm(successful_confs, desc="pMHC conf", position=0, disable=disable_progress_bar):

			numTries = 1
			best_energy = float("inf")
			pMHC_complex = pMHC(pdb_filename=filestore + "/5_openMM_conformations/13_connected_pMHC_complexes/pMHC_" + conf_index + ".pdb", 
								peptide=peptide)
			for minimization_effort in tqdm(range(1, numTries + 1),  desc="No. of tries", position=1,
											leave=leave_progress_bar, disable=disable_progress_bar):
				best_energy = pMHC_complex.minimizeConf(filestore, best_energy, no_constraints_openmm, device)
				with open(filestore + "/5_openMM_conformations/05_per_peptide_results/peptide_" + conf_index + ".log", 'w') as peptide_handler:
					peptide_handler.write(conf_index + ",Successfully Modeled," + str(best_energy) + "\n")

		# Rescoring and re-filtering resulting conformations
		if verbose: print("\nRescoring and re-filtering resulting conformations:")
		#arg_list = list(map(lambda conf_index: (conf_index, filestore, rcd_num_loops, peptide_template_anchors_xyz, anchor_tol, tolerance_anchors, min_with_smina), successful_confs))
		arg_list = []
		for conf_index in successful_confs:
			arg_list.append((conf_index, filestore, rcd_num_loops, anchor_filtering_data_dict[int(conf_index[0])][0], 
							 anchor_tol, anchor_filtering_data_dict[int(conf_index[0])][1], min_with_smina))
		with WorkerPool(n_jobs=min(num_cores, len(successful_confs))) as pool:
			results = pool.map(rescoring_after_openmm, arg_list, progress_bar=verbose)

		copy_batch_of_files(filestore + '/5_openMM_conformations/10_pMHC_complexes/',
							filestore + '/6_final_conformations/',
							query="pMHC_")

		best_conf_dir = filestore + '/5_openMM_conformations'
		if verbose: print("\n\nEnd of OpenMM step !!!")
		create_csv_from_list_of_files(filestore + '/5_openMM_conformations/total_results.csv', glob.glob(filestore + '/5_openMM_conformations/05_per_peptide_results/*.log'))
		results_csv = pretty_print_analytics(filestore + '/5_openMM_conformations/total_results.csv', verbose=verbose)
		results_csv.to_csv(filestore + '/5_openMM_conformations/successful_conformations_statistics.csv', index=False)
		results_csv.to_csv(temp_files_storage + '/successful_conformations_statistics.csv', index=False)	
	else:
		initialize_dir(filestore + '/5_final_conformations/')
		copy_batch_of_files(filestore + '/4_SMINA_data/10_pMHC_complexes/',
							filestore + '/5_final_conformations',
							query="pMHC_")

	# Control whether there are no conformations. If they do, store the best one and continue.
	# If not, either abort or force restart (for round one)
	if(results_csv.shape[0] == 0):
		print('No conformations were produced...')
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

	# Delete intermediate files if flag is true
	if not keep_all_files:
		dir_list = ['/4_SMINA_data', '/3_RCD_data', '/2_input_to_RCD', '/1_alignment_files']
		if(score_with_openmm and results_csv.shape[0] > 0):
			dir_list.append('/5_openMM_conformations')
		remove_dirs([filestore + dir for dir in dir_list])

	print("\n\nEnd of APE-Gen")

if __name__ == "__main__":
	apegen(sys.argv[1:])
