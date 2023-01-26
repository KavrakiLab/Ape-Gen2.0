import os
from sys import stdout
import shutil
import sys
import re
from pathlib import Path

import pandas as pd
import numpy as np

from pdbtools import pdb_merge, pdb_tidy, pdb_reatom, pdb_sort, pdb_selchain, pdb_rplchain, pdb_selmodel
from Bio.Align import substitution_matrices

from constraint import *

from nltk import ngrams

# PDBFIXER
from pdbfixer import PDBFixer
from openmm.app import PDBFile

## MACROS

# set verbose
global verbose_var
verbose_var = False
def verbose():
	return verbose_var
def set_verbose(val):
	global verbose_var
	verbose_var = val

# Three-to-one (and vice versa) AA transformation
standard_three_to_one_letter_code = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', \
									 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', \
									 'GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', 'ILE':'I', \
									 'LEU':'L', 'MET':'M', 'PHE':'F', 'TYR':'Y', 'TRP':'W'}
standard_one_to_three_letter_code = {value: key for key, value in standard_three_to_one_letter_code.items()}


non_standard_three_to_one_letter_code = {'SEP':'pS', 'TPO':'pT', 'PTR':'pY', 'mmK':'m', 'mdK':'d', 'mtK':'t', 'xhC':'h', 'xoC':'o', 'xdC':'d'}
non_standard_one_to_three_letter_code = {value: key for key, value in non_standard_three_to_one_letter_code.items()}

all_three_to_one_letter_codes = dict(standard_three_to_one_letter_code, **non_standard_three_to_one_letter_code)
all_one_to_three_letter_codes = dict(standard_one_to_three_letter_code, **non_standard_one_to_three_letter_code)


anchor_dictionary = {'8': {'1': 'N_+0', '2': 'N_+1', '3': 'N_+2', '4': 'N_+3',
						   '5': 'C_-3', '6': 'C_-2', '7': 'C_-1', '8': 'C_+0'},
					 '9': {'1': 'N_+0', '2': 'N_+1', '3': 'N_+2', '4': 'N_+3', '5': 'N_+4',
						   '6': 'C_-3', '7': 'C_-2', '8': 'C_-1', '9': 'C_+0'},
					 '10': {'1': 'N_+0', '2': 'N_+1', '3': 'N_+2', '4': 'N_+3', '5': 'N_+4',
						   '6': 'C_-4', '7': 'C_-3', '8': 'C_-2', '9': 'C_-1', '10': 'C_+0'},
					 '11': {'1': 'N_+0', '2': 'N_+1', '3': 'N_+2', '4': 'N_+3', '5': 'N_+4', '6': 'N_+5',
						   '7': 'C_-4', '8': 'C_-3', '9': 'C_-2', '10': 'C_-1', '11': 'C_+0'},
					 '12': {'1': 'N_+0', '2': 'N_+1', '3': 'N_+2', '4': 'N_+3', '5': 'N_+4', '6': 'N_+5',
						   '7': 'C_-5', '8': 'C_-4', '9': 'C_-3', '10': 'C_-2', '11': 'C_-1', '12': 'C_+0'},
					 '13': {'1': 'N_+0', '2': 'N_+1', '3': 'N_+2', '4': 'N_+3', '5': 'N_+4', '6': 'N_+5', '7': 'N_+6',
						   '8': 'C_-5', '9': 'C_-4', '10': 'C_-3', '11': 'C_-2', '12': 'C_-1', '13': 'C_+0'},
					 '14': {'1': 'N_+0', '2': 'N_+1', '3': 'N_+2', '4': 'N_+3', '5': 'N_+4', '6': 'N_+5', '7': 'N_+6',
						   '8': 'C_-6', '9': 'C_-5', '10': 'C_-4', '11': 'C_-3', '12': 'C_-2', '13': 'C_-1', '14': 'C_+0'},
					 '15': {'1': 'N_+0', '2': 'N_+1', '3': 'N_+2', '4': 'N_+3', '5': 'N_+4', '6': 'N_+5', '7': 'N_+6', '8': 'N_+7',
						   '9': 'C_-6', '10': 'C_-5', '11': 'C_-4', '12': 'C_-3', '13': 'C_-2', '14': 'C_-1', '15': 'C_+0'}}

rev_anchor_dictionary = {'N_+0' : {'8': 1, '9': 1, '10': 1, '11': 1, '12': 1, '13': 1, '14': 1, '15': 1},
						 'N_+1' : {'8': 2, '9': 2, '10': 2, '11': 2, '12': 2, '13': 2, '14': 2, '15': 2},
						 'N_+2' : {'8': 3, '9': 3, '10': 3, '11': 3, '12': 3, '13': 3, '14': 3, '15': 3},
						 'N_+3' : {'8': 4, '9': 4, '10': 4, '11': 4, '12': 4, '13': 4, '14': 4, '15': 4},
						 'N_+4' : {'8': 5, '9': 5, '10': 5, '11': 5, '12': 5, '13': 5, '14': 5, '15': 5},
						 'N_+5' : {'8': 6, '9': 6, '10': 6, '11': 6, '12': 6, '13': 6, '14': 6, '15': 6},
						 'N_+6' : {'8': 7, '9': 7, '10': 7, '11': 7, '12': 7, '13': 7, '14': 7, '15': 7},
						 'N_+7' : {'8': 8, '9': 8, '10': 8, '11': 8, '12': 8, '13': 8, '14': 8, '15': 8},
						 'C_-6' : {'8': 2, '9': 3, '10': 4, '11': 5, '12': 6, '13': 7, '14': 8, '15': 9},
						 'C_-5' : {'8': 3, '9': 4, '10': 5, '11': 6, '12': 7, '13': 8, '14': 9, '15': 10},
						 'C_-4' : {'8': 4, '9': 5, '10': 6, '11': 7, '12': 8, '13': 9, '14': 10, '15': 11},
						 'C_-3' : {'8': 5, '9': 6, '10': 7, '11': 8, '12': 9, '13': 10, '14': 11, '15': 12},
						 'C_-2' : {'8': 6, '9': 7, '10': 8, '11': 9, '12': 10, '13': 11, '14': 12, '15': 13},
						 'C_-1' : {'8': 7, '9': 8, '10': 9, '11': 10, '12': 11, '13': 12, '14': 13, '15': 14},
						 'C_+0' : {'8': 8, '9': 9, '10': 10, '11': 11, '12': 12, '13': 13, '14': 14, '15': 15}}

## FUNCTIONS

def initialize_dir(dir_name):
	if type(dir_name) == str:
		if os.path.exists(dir_name):
			for root, dirs, files in os.walk(dir_name):
				for f in files:
					os.unlink(os.path.join(root, f))
				for d in dirs:
					shutil.rmtree(os.path.join(root, d))
		else:
			os.umask(0)
			os.makedirs(dir_name,mode=0o777)
	elif type(dir_name) == list:
		for dir in dir_name:
			initialize_dir(dir)
		
def copy_file(src, dst):
	shutil.copy(src, dst)

def move_file(src, dst):
	shutil.move(src, dst)

def move_batch_of_files(src, dst, query):
	files = os.listdir(src)
	for f in files:
		if (query in f): shutil.move(src + f, dst)

def copy_batch_of_files(src, dst, query):
	files = os.listdir(src)
	for f in files:
		if (query in f): shutil.copy(src + f, dst)

def remove_file(filename):
	os.remove(filename)
	
def add_sidechains(pdb_filename, filestore, add_hydrogens="Yes", peptide_idx=-1, remove_heterogens=True,
				   add_solvent=False, keep_IDs=False):

	fixer = PDBFixer(filename=pdb_filename)
	fixer.findMissingResidues()
	if remove_heterogens: fixer.removeHeterogens(True) #  True keeps water molecules while removing all other heterogens, REVISIT!
	fixer.findMissingAtoms()
	fixer.addMissingAtoms()
	if add_hydrogens != "none": fixer.addMissingHydrogens(7.0)
	if add_solvent: fixer.addSolvent(fixer.topology.getUnitCellDimensions())
	if peptide_idx != -1:
			pdb_filename = filestore + '/02_add_sidechains/PTMed_' + str(peptide_idx) + '.pdb'
	
	fp = open(pdb_filename, 'w')
	PDBFile.writeFile(fixer.topology, fixer.positions, fp, keepIds=keep_IDs)
	fp.close()

	if peptide_idx != -1:
		copy_file(filestore + '/02_add_sidechains/PTMed_' + str(peptide_idx) + '.pdb', 
							filestore + '/03_PTMed_peptides/PTMed_' + str(peptide_idx) + '.pdb')
	return pdb_filename

def add_missing_residues(pdb_filename, filestore):
	fixer = PDBFixer(filename=pdb_filename)
	fixer.findMissingResidues()
	fixer.findMissingAtoms()
	fixer.addMissingAtoms()
	fp = open(pdb_filename, 'w')
	PDBFile.writeFile(fixer.topology, fixer.positions, fp, keepIds=True)
	fp.close()
	return pdb_filename

def apply_mutations(pdb_filename, filestore, mutation_list):
	fixer = PDBFixer(filename=pdb_filename)
	fixer.applyMutations(mutation_list, "C")
	fixer.findMissingResidues()
	fixer.findMissingAtoms()
	fixer.addMissingAtoms()
	fp = open(pdb_filename, 'w')
	PDBFile.writeFile(fixer.topology, fixer.positions, fp, keepIds=True)
	fp.close()
	return pdb_filename

def merge_and_tidy_pdb(list_of_pdbs, dst):
	merged = pdb_merge.run(pdb_merge.check_input(list_of_pdbs))
	sorteded = pdb_sort.run(merged, sorting_keys='-RC')
	tidied = pdb_tidy.run(sorteded, strict=True)
	reatomed = pdb_reatom.run(tidied, starting_value=1)
	with open(dst, 'w') as pdb_file:
		pdb_file.write(''.join(reatomed))
	pdb_file.close()

def apply_function_to_file(func, input_filename, output_filename="", **kwargs):
	if output_filename == "": output_filename = input_filename

	overwritten = func(input_filename, **{key: value for key, value in kwargs.items() if key in func.__code__.co_varnames})
	overwritten = ''.join(overwritten)

	with open(output_filename, 'w') as output:
		output.write(overwritten)
	return output_filename

def rename_chains(pdb_filename, chain_from, chain_to, dst):
	renamed = pdb_rplchain.run(open(pdb_filename, 'r'), (chain_from, chain_to))
	with open(dst, 'w') as file:
		file.write(''.join(renamed))
	file.close()

def select_models(pdb_filename, best_indexes, dst):
	selected = pdb_selmodel.run(open(pdb_filename, 'r'), model_set=best_indexes)
	with open(dst, 'w') as file:
		file.write(''.join(selected))
	file.close()

def split_receptor_and_peptide(pdb_file):
	pdb_path, pdb_filename = os.path.split(pdb_file)
	pdb_filename = pdb_filename.replace(".pdb", "")
	receptor_filename = pdb_path + '/' + pdb_filename + "_receptor.pdb"
	peptide_filename = pdb_path + '/' + pdb_filename + "_peptide.pdb"

	receptor = pdb_selchain.run(open(pdb_file, 'r'), ('A', 'B'))
	with open(receptor_filename, 'w') as file:
		file.write(''.join(receptor))
	file.close()

	peptide = pdb_selchain.run(open(pdb_file, 'r'), ('C',))
	with open(peptide_filename, 'w') as file:
		file.write(''.join(peptide))
	file.close()

	return (receptor_filename, peptide_filename)

def filter_chains(pdb_file, chains, dst):
	filtered = pdb_selchain.run(open(pdb_file, 'r'), chains)
	with open(dst, 'w') as file:
		file.write(''.join(filtered))
	file.close()

def remove_remarks_and_others_from_pdb(pdb_file, records=('ATOM', 'TER', 'END ')): 
	fhandle = open(pdb_file, 'r')
	for line in fhandle:
		if line.startswith(records):
			yield line
	fhandle.close()

def replace_chains(pdb_file, chain_from, chain_to):
	records = ('ATOM', 'HETATM', 'ANISOU')
	fhandle = open(pdb_file, 'r')
	for line in fhandle:
		if line.startswith(records):
			if line[21] == chain_from:
				yield line[:21] + chain_to + line[22:]
				continue
		yield line
	fhandle.close()

def replace_HETATM(pdb_file):
	fhandle = open(pdb_file, 'r')
	for line in fhandle:
		if line.startswith("HETATM"):
			yield "ATOM  " + line[6:]
			continue
		yield line
	fhandle.close()

def delete_elements(pdb_file, element_set, chains):
	records = ('ATOM', 'HETATM', 'ANISOU')
	fhandle = open(pdb_file, 'r')
	for line in fhandle:
		if line.startswith(records):
			if line[13:15].strip() in element_set and line[21].strip() in chains:
				continue
		yield line
	fhandle.close()

def count_number_of_atoms(pdb_file):
	fhandle = open(pdb_file, 'r')
	count = 0
	for line in fhandle:
		if line.startswith('ATOM'):
			count +=1
	fhandle.close()
	return count

def create_csv_from_list_of_files(csv_filename, list_of_files):
	with open(csv_filename, 'wb') as outfile:
		for filename in list_of_files:
			with open(filename, 'rb') as readfile:
				shutil.copyfileobj(readfile, outfile)

def pretty_print_analytics(csv_location, verbose=True):

	results_csv = pd.read_csv(csv_location, names=['Round', 'Peptide index', 'Debug', 'Affinity'])
	results_csv.sort_values(by=['Round', 'Peptide index'], inplace=True)
	results_csv.to_csv(csv_location, index=False)

	results_csv = results_csv[results_csv['Affinity'] != '-']
	no_of_final_conformations = results_csv.shape[0]
	if verbose:
		print("Total number of overall conformations: " + str(no_of_final_conformations))
		print("\nAnalytics:")
		print(results_csv.to_markdown(index=False))
	return results_csv

def score_sequences(seq1, seq2, matrix, gap_penalty, norm):
	score = 0
	num_gaps = 0
	for A,B in zip(seq1,seq2):
		diag = ('-'==A) or ('-'==B)
		score += 0 if diag else matrix[A,B]
		num_gaps += 1 if diag else 0
	score = score/norm
	return score + num_gaps*gap_penalty

def anchor_alignment(seq1, seq2, anchor_diff1, anchor_diff2):
	extra_in_beginning = ''.join('-'*abs(anchor_diff1))
	extra_in_end = ''.join('-'*abs(anchor_diff2))
	if anchor_diff1 > 0:
		temp_seq1 = extra_in_beginning + seq1
		temp_seq2 = seq2
	else:
		temp_seq1 = seq1
		temp_seq2 = extra_in_beginning + seq2
	if anchor_diff2 > 0:
		temp_seq1 = temp_seq1 + extra_in_end
	else:
		temp_seq2 = temp_seq2 + extra_in_end
	return temp_seq1, temp_seq2

def calculate_anchors_given_alignment(seq, seq_template, anchor_1, anchor_2):

	extra_in_beginning = len(seq)-len(seq.lstrip('-'))
	extra_in_beginning_template = len(seq_template)-len(seq_template.lstrip('-'))
	extra_in_end = len(seq)-len(seq.rstrip('-'))
	extra_in_end_template = len(seq_template)-len(seq_template.rstrip('-'))
	anchor_1 = anchor_1 - extra_in_beginning + extra_in_beginning_template
	anchor_2 = anchor_2 - extra_in_beginning + extra_in_beginning_template - extra_in_end + extra_in_end_template 
	return (anchor_1, anchor_2)

## RESIDUE RENAMING AFTER SMINA FLEXIBILITY OUTPUT

def csp_solver(edge_list, residue, atom_indexes, CA_loc, C_loc):

	# Note to change the 4-letter atoms if need be!
	atom_dict = {'ALA':[["CA", "C", "CB"]],
				 'VAL':[["CA", "C", "CB", "CG1", "CG2"]],
				 'ILE':[["CA", "C", "CB", "CG1", "CG2", "CD"]],
				 'LEU':[["CA", "C", "CB", "CG", "CD1", "CD2"]],
				 'MET':[["CA", "C", "CB", "CG", "SD", "CE"]],
				 'PHE':[["CA", "C", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
				 'TYR':[["CA", "C", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "HH"]],
				 'TRP':[["CA", "C", "CB", "CG", "CD1", "CD2", "NE1", "HE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]],
				 'SER':[["CA", "C", "CB", "OG", "HG1"]],
				 'THR':[["CA", "C", "CB", "OG1", "HG1", "CG2"]],
				 'ASN':[["CA", "C", "CB", "CG", "OD1", "ND2", "HD21", "HD22"]],
				 'GLN':[["CA", "C", "CB", "CG", "CD", "OE1", "NE2", "HE21", "HE22"]],
				 'CYS':[["CA", "C", "CB", "SG"]],
				 'GLY':[["CA", "C"]],
				 'PRO':[["CA", "C", "CB", "CG", "CD"]],
				 'ARG':[["CA", "C", "CB", "CG", "CD", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"]],
				 'HIS':[["CA", "C", "CB", "CG", "ND1", "HD1", "CE1", "NE2", "CD2"], 
						["CA", "C", "CB", "CG", "ND1", "HE2", "CE1", "NE2", "CD2"]],
				 'LYS':[["CA", "C", "CB", "CG", "CD", "CE", "NZ", "HZ1", "HZ2", "HZ3"]],
				 'ASP':[["CA", "C", "CB", "CG", "OD1", "OD2"]],
				 'GLU':[["CA", "C", "CB", "CG", "CD", "OE1", "OE2"]]
				}

	constraint_dict = {'ALA':[[["CA", "C"], ["CA", "CB"]]],
					   'VAL':[[["CA", "C"], ["CA", "CB"], ["CB", "CG1"], ["CB", "CG2"]]],
					   'ILE':[[["CA", "C"], ["CA", "CB"], ["CB", "CG1"], ["CB", "CG2"], ["CG1", "CD"]]],
					   'LEU':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD1"], ["CG", "CD2"]]],
					   'MET':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "SD"], ["SD", "CE"]]],
					   'PHE':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD1"], ["CG", "CD2"], ["CD1", "CE1"], ["CD2", "CE2"], ["CE1", "CZ"], ["CE2", "CZ"]]],
					   'TYR':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD1"], ["CG", "CD2"], ["CD1", "CE1"], ["CD2", "CE2"], ["CE1", "CZ"], ["CE2", "CZ"], ["CZ", "OH"], ["OH", "HH"]]],
					   'TRP':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD1"], ["CG", "CD2"], ["CD1", "NE1"], ["NE1", "CE2"], ["NE1", "HE1"], ["CD2", "CE2"], ["CD2", "CE3"], ["CE2", "CZ2"], ["CE3", "CZ3"], ["CZ2", "CH2"], ["CH2", "CZ3"]]],
					   'SER':[[["CA", "C"], ["CA", "CB"], ["CB", "OG"], ["OG", "HG1"]]],
					   'THR':[[["CA", "C"], ["CA", "CB"], ["CB", "OG1"], ["OG1", "HG1"], ["CB", "CG2"]]],
					   'ASN':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "ND2"], ["ND2", "HD21"], ["ND2", "HD22"]]],
					   'GLN':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "OE1"], ["CD", "NE2"], ["NE2", "HE21"], ["NE2", "HE22"]]],
					   'CYS':[[["CA", "C"], ["CA", "CB"], ["CB", "SG"]]],
					   'GLY':[[["CA", "C"]]],
					   'PRO':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"]]],
					   'ARG':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "NE"], ["NE", "HE"], ["NE", "CZ"], ["CZ", "NH1"], ["CZ", "NH2"], ["NH1", "HH11"], ["NH1", "HH12"], ["NH2", "HH21"], ["NH2", "HH22"]]],
					   'HIS':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "ND1"], ["ND1", "CE1"], ["ND1", "HD1"], ["CE1", "NE2"], ["NE2", "CD2"], ["CD2", "CG"]],
							  [["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "ND1"], ["ND1", "CE1"], ["NE2", "HE2"], ["CE1", "NE2"], ["NE2", "CD2"], ["CD2", "CG"]]],
					   'LYS':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "CE"], ["CE", "NZ"], ["NZ", "HZ1"], ["NZ", "HZ2"], ["NZ", "HZ3"]]],
					   'ASP':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "OD1"], ["CG", "OD2"]]],
					   'GLU':[[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "OE1"], ["CD", "OE2"]]]
				}

	no_of_cases = len(atom_dict[residue])
	for i in range(no_of_cases):
		# Formulating the atom matching problem as a CSP:
		problem = Problem()

		# Adding atoms as variables to the CSP (along with possible index values they can take (CA and C are known):
		problem.addVariable("CA", [CA_loc])
		problem.addVariable("C", [C_loc])
		problem.addVariables((atom_dict[residue][i])[2:], atom_indexes)

		# Adding Constraints to the CSP
		# 1. No atoms share the same index
		problem.addConstraint(AllDifferentConstraint())

		# 2. Topological constraints that are associated with each residue:
		for elem in constraint_dict[residue][i]:
			problem.addConstraint(lambda a, b: ([a,b] in edge_list) and ([b,a] in edge_list), elem)

		# 3. Find problem solution
		if len(problem.getSolutions()) > 0:
			solution = problem.getSolutions()[0]
			return pd.DataFrame(data={'atom_name': solution.keys(), 'atom_number': list(solution.values())})

	return pd.DataFrame(data={})

def extract_CONECT_from_pdb(pdb_file):

	edge_list = []
	taken = remove_remarks_and_others_from_pdb(pdb_file, records=("CONECT "))
	conect_fields = ''.join(taken)
	for line in conect_fields.split('\n'):
		cleaned_line = [elem for elem in line.strip().split(' ') if len(elem) > 0]
		try:
			pivot = cleaned_line[1]
			for atom_index in cleaned_line[2:]:
				edge_list.append([int(pivot), int(atom_index)])
		except IndexError:
			pass
	return edge_list

## PREPARE PTMs FOR OPENMM

def replace_CONECT_fields(pdb_file, index_df, external_bonds_list):

	# Take care of the external bonds first
	for bond in external_bonds_list:
		yield "CONECT " + str(bond[0]) + " " + str(bond[1]) + "\n"

	# Now the rest of the bonds
	taken = remove_remarks_and_others_from_pdb(pdb_file, records=("CONECT "))
	conect_fields = ''.join(taken)
	for line in conect_fields.split('\n'):
		if line == "":
			yield ""
		else:
			cleaned_line = [str(index_df[index_df['atom_name'] == elem]['atom_number'].item()) for elem in line.strip().split(' ')[1:] if len(elem) > 0]
			cleaned_line.append('\n')
			yield "CONECT " + ' '.join(cleaned_line)

def merge_connect_fields(list_of_pdbs, dst):
	merged = pdb_merge.run(pdb_merge.check_input(list_of_pdbs))
	sorteded = pdb_sort.run(merged, sorting_keys='-RC') # Potentially breaking -> Not working?
	with open(dst, 'w') as pdb_file:
		pdb_file.write(''.join(sorteded))
	pdb_file.close()

def extract_anchors(peptide, MHC, frequencies):

	# Preprocessing for feature extration
	frequencies = frequencies[frequencies['length'].isin([8,9,10,11])]
	frequencies_first_part = frequencies[frequencies['position'].isin([1,2,3])]
	frequencies_second_part = frequencies[frequencies['length'] - frequencies['position'] <= 2].copy()
	frequencies_second_part['position'] = frequencies_second_part['position'].sub(frequencies_second_part['length'], axis = 0)
	frequencies_first_part = frequencies_first_part.groupby(['allele', 'position']).max().reset_index().drop(['length', 'cutoff_fraction', 'cutoff_count'], axis=1)
	frequencies_second_part = frequencies_second_part.groupby(['allele', 'position']).max().reset_index().drop(['length', 'cutoff_fraction', 'cutoff_count'], axis=1)

	# Feature for pos 3 only
	aa_volume = {'G': 60.1, 'A': 88.6, 'S': 89.0, 'C': 108.5, 'D' : 111.1, 'P' : 112.7, 'N': 114.1, 'T': 116.1,
				 'E': 138.4, 'V': 140.0, 'Q': 143.8, 'H': 153.2, 'M': 162.9, 'I': 166.7, 'L': 116.7, 'K': 168.6,
				 'R': 173.4, 'F': 189.9, 'Y': 193.6, 'W': 227.8}

	first_part_of_peptide = peptide[:3]
	pep_sequence = list(first_part_of_peptide)
	freq_features = frequencies_first_part[(frequencies_first_part['allele'] == MHC)]
	potential_pos_1 = freq_features[freq_features['position'] == 2][pep_sequence[0]].values[0]
	potential_pos_3 = freq_features[freq_features['position'] == 2][pep_sequence[2]].values[0]
	stability_pos_2 = freq_features[freq_features['position'] == 2][pep_sequence[1]].values[0]
	inertia_pos_1 = freq_features[freq_features['position'] == 1][pep_sequence[0]].values[0]
	inertia_pos_3 = freq_features[freq_features['position'] == 3][pep_sequence[0]].values[0]
	volume_pos_3 = aa_volume[pep_sequence[0]] + aa_volume[pep_sequence[1]]

	second_part_of_peptide = peptide[-3:]
	pep_sequence = list(second_part_of_peptide)
	freq_features = frequencies_second_part[(frequencies_second_part['allele'] == MHC)]
	potential_pos_C2 = freq_features[freq_features['position'] == 0][pep_sequence[0]].values[0]
	potential_pos_C1 = freq_features[freq_features['position'] == 0][pep_sequence[1]].values[0]
	stability_pos_C = freq_features[freq_features['position'] == 0][pep_sequence[2]].values[0]
	inertia_pos_C2 = freq_features[freq_features['position'] == -2][pep_sequence[0]].values[0]
	inertia_pos_C1 = freq_features[freq_features['position'] == -1][pep_sequence[1]].values[0]

	# N-termini
	anchor_1 = "2"
	if (potential_pos_3 > 0.195) and (volume_pos_3 < 180) and (stability_pos_2 < 0.1):
		anchor_1 = "3"
	if (potential_pos_1 > 0.195) and (inertia_pos_1 < 0.14) and (stability_pos_2 < 0.08):
		anchor_1 = "1"

	# C-termini
	anchor_2 = str(len(peptide))
	if (potential_pos_C1 > 0.16) and (inertia_pos_C1 < 0.14) and (stability_pos_C < 0.08):
		anchor_2 = str(len(peptide) - 1) 
	if (potential_pos_C2 > 0.25) and (stability_pos_C < 0.02):
		anchor_2 = str(len(peptide) - 2) 

	return ",".join([anchor_1, anchor_2])

def extract_anchors_PMBEC(peptide, MHC, frequencies):
	
	# Load and set the SMM matrix
	smm_matrix = pd.read_csv('./helper_files/PMBEC/' + MHC + '.txt', sep = '\t', skiprows=1, header = None, nrows=20).transpose()
	new_header = smm_matrix.iloc[0] #grab the first row for the header
	smm_matrix = smm_matrix[1:] #take the data less the header row
	smm_matrix.columns = new_header #set the header row as the df header
	smm_matrix = smm_matrix.reset_index().rename(columns={'index': 'position'})

	matrix = frequencies[(frequencies['allele'] == MHC) & (frequencies['length'] == 9)]
	
	# N-termini_candidate_1
	first_part_of_peptide = peptide[:3]
	pep_sequence = list(first_part_of_peptide)

	potential_pos_11 = smm_matrix[smm_matrix['position'] == 2][pep_sequence[0]].values[0]
	inertia_pos_11 = smm_matrix[smm_matrix['position'] == 1][pep_sequence[0]].values[0]
	will_11 = potential_pos_11 - inertia_pos_11
	potential_12 = smm_matrix[smm_matrix['position'] == 3][pep_sequence[1]].values[0]
	inertia_12 = smm_matrix[smm_matrix['position'] == 2][pep_sequence[1]].values[0]
	will_12 = potential_12 - inertia_12
	potential_pos_13 = smm_matrix[smm_matrix['position'] == 4][pep_sequence[2]].values[0]
	inertia_pos_13 = smm_matrix[smm_matrix['position'] == 3][pep_sequence[2]].values[0]
	will_13 = potential_pos_13 - inertia_pos_13
	
	# N-termini_candidate_3
	first_part = matrix[(matrix['position'].isin([1,2,3,4,5]))]
	anchor_1 = first_part.iloc[:, 5:].max(1).argmax() + 1
	if anchor_1 in [1, 4, 5]: # This is because it probably cannot happen
		anchor_1 = 2
	
	first_part_of_peptide = peptide[:(anchor_1 + 2)]
	pep_sequence = list(first_part_of_peptide)

	inertia_pos_31 = smm_matrix[smm_matrix['position'] == 1][pep_sequence[0]].values[0]
	potential_32 = smm_matrix[smm_matrix['position'] == anchor_1 - 1][pep_sequence[anchor_1 - 1]].values[0]
	inertia_32 = smm_matrix[smm_matrix['position'] == anchor_1][pep_sequence[anchor_1 - 1]].values[0]
	will_32 = potential_32 - inertia_32
	potential_pos_33 = smm_matrix[smm_matrix['position'] == anchor_1][pep_sequence[anchor_1]].values[0]
	inertia_pos_33 = smm_matrix[smm_matrix['position'] == anchor_1 + 1][pep_sequence[anchor_1]].values[0]
	will_33 = potential_pos_33 - inertia_pos_33
	
	# C-termini_candidate
	second_part_of_peptide = peptide[7:]
	second_part = matrix[(matrix['position'] >= 8) & (matrix['position'] <= len(peptide))]
	pep_sequence = list(second_part_of_peptide)
	stability_c = smm_matrix[smm_matrix['position'] == 9][pep_sequence[len(peptide) - 8]].values[0]
	min_will_c = float("inf")
	arg_will_c = 0
	for pos in range(8, len(peptide)):
		potential_pos_c = smm_matrix[smm_matrix['position'] == 9][pep_sequence[pos - 8]].values[0]
		intertia_pos_c = smm_matrix[smm_matrix['position'] == 8][pep_sequence[pos - 8]].values[0]
		will_pos_c = potential_pos_c - intertia_pos_c
		if min_will_c > will_pos_c:
			min_will_c = will_pos_c
			arg_will_c = pos
	print(stability_c)
	print(min_will_c)
	anchor_1 = "2"
	if (will_11 < -0.3) and (will_12 < 0.05) and (will_13 < -0.06):
		anchor_1 = "1"
	if (inertia_pos_31 > 0.0) and (will_32 < 0.75) and (will_33 < -0.35):
		anchor_1 = "3"
	anchor_2 = str(len(peptide))
	if (stability_c > 0.0) and (min_will_c < -0.5):
		anchor_2 = str(arg_will_c)
	
	return ",".join([anchor_1, anchor_2])

def predict_anchors(peptide, MHC):

	frequencies = pd.read_csv("./helper_files/mhcflurry.ba.frequency_matrices.csv")
	frequencies = frequencies[(frequencies['cutoff_fraction'] == 0.01)]
	frequencies['X'] = np.zeros(frequencies.shape[0])
	frequencies_alleles = pd.unique(frequencies['allele'])

	if MHC in frequencies_alleles:
		anchors = extract_anchors(peptide, MHC, frequencies)
		anchor_status = "Known"
	else:
		anchors = "2," + str(len(peptide))
		anchor_status = "Not Known"
	return anchors, anchor_status

def predict_anchors_PMBEC(peptide, MHC):

	frequencies = pd.read_csv("./helper_files/mhcflurry.ba.frequency_matrices.csv")
	frequencies = frequencies[(frequencies['cutoff_fraction'] == 0.01)]
	frequencies['X'] = np.zeros(frequencies.shape[0])

	frequencies_alleles = os.listdir('./helper_files/PMBEC/') 
	frequencies_alleles = [x.split('.')[0] for x in frequencies_alleles]
	
	if (len(peptide) <= 7) or MHC not in frequencies_alleles:
		anchors = "2," + str(len(peptide))
		anchor_status = "Not Known"
	else:
		anchors = extract_anchors_PMBEC(peptide, MHC, frequencies)
		anchor_status = "Known"
	return anchors, anchor_status
