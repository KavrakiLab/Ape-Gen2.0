import os
import shutil
import sys
import re
from pathlib import Path

import pandas as pd
import numpy as np

from pdbtools import pdb_merge, pdb_tidy, pdb_reatom, pdb_sort, pdb_selchain, pdb_rplchain
from Bio.Align import substitution_matrices

from constraint import *

from nltk import ngrams

## MACROS

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

def AA_error_checking(amino_acid):
	if (amino_acid not in standard_three_to_one_letter_code.values()) and (amino_acid not in non_standard_three_to_one_letter_code.values()):
		print("The provided amino acid in the sequence is wrong")
		sys.exit(0)

def process_anchors(anchors, pep_seq):
	# Returns the set of anchors
	pep_length = len(re.sub('[a-z]', '', pep_seq)) # Remove any PTMs that may still exist in the sequence
	anchor_not = set([anchor_dictionary[str(pep_length)][str(aa_index)] for aa_index in anchors.split(",")])
	return anchor_not

def jaccard_distance(a, b):
	# Computes jaccard distance between 2 sets
	c = a.intersection(b)
	return float(len(c)) / (len(a) + len(b) - len(c))

def initialize_dir(*dir_names):
	for dir_name in dir_names:
		if os.path.exists(dir_name):
			shutil.rmtree(dir_name)
		else:
			pass
		os.makedirs(dir_name)

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

def merge_and_tidy_pdb(list_of_pdbs, dst):
	merged = pdb_merge.run(pdb_merge.check_input(list_of_pdbs))
	sorteded = pdb_sort.run(merged, sorting_keys='-RC') # Potentially breaking -> Not working?
	tidied = pdb_tidy.run(sorteded, strict=True)
	reatomed = pdb_reatom.run(tidied, starting_value=1)
	with open(dst, 'w') as pdb_file:
		pdb_file.write(''.join(reatomed))
	pdb_file.close()

def apply_function_to_file(func, filename, **kwargs):
	overwritten = func(filename, **{key: value for key, value in kwargs.items() if key in func.__code__.co_varnames})
	overwritten = ''.join(overwritten)

	with open(filename, 'w') as altered:
		altered.write(overwritten)

def rename_chains(pdb_filename, chain_from, chain_to, dst):
	renamed = pdb_rplchain.run(open(pdb_filename, 'r'), (chain_from, chain_to))
	with open(dst, 'w') as file:
		file.write(''.join(renamed))
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
	   
def create_csv_from_list_of_files(csv_filename, list_of_files):
	with open(csv_filename, 'wb') as outfile:
		for filename in list_of_files:
			with open(filename, 'rb') as readfile:
				shutil.copyfileobj(readfile, outfile)

def pretty_print_analytics(csv_location):

	results_csv = pd.read_csv(csv_location, names=['Round', 'Peptide index', 'Debug', 'Affinity'])
	results_csv.sort_values(by=['Round', 'Peptide index'], inplace = True)
	results_csv.to_csv(csv_location, index=False)

	results_csv = results_csv[results_csv['Affinity'] != '-']
	no_of_final_conformations = results_csv.shape[0]
	print("Total number of overall conformations: " + str(no_of_final_conformations))
	print("\nAnalytics:")
	print(results_csv.to_markdown(index = False))
	return results_csv

## PTMs

# Different PTMs

phosphorylation_list = ['pS', 'pT', 'pY']
acetylation_list = ['aK'] # Check details on pytms -> This is nmot working well, it renames all hydrogens to PyMOL ones
carbamylation_list = ['cK'] # Check details on pytms
citrullination_list = ['cR'] 
methylation_list = ['mmK', 'mdK', 'mtK'] # Check details on pytms
nitration_list = ['nY', 'nW'] # Check details on pytms
s_nitrosylation_list = ['nC']
p_hydroxylation_list = ['nP'] # Check details on pytms
# malondialdehyde_list = ['maK'] # Check details on pytms, but probably it's too complicated to use that one
c_oxidation_list = ['xhC', 'xoC', 'xdC'] # Check details on pytms
m_oxidation_list = ['oM'] # Check details on pytms


def sequence_PTM_processing(sequence):
	sequence_list = re.sub( r"([A-Z])", r"\1 ", sequence).split() # Split pep sequence while retaining the PTM
		
	PTM_list = []
	for i, amino_acid in enumerate(sequence_list):
		if len(amino_acid) > 1:
			prefix = PTM_error_checking(amino_acid)
			AA_error_checking(amino_acid[1])
			PTM_list.append(prefix + str(i + 1))
		else:
			AA_error_checking(amino_acid)
	return PTM_list

def PTM_error_checking(amino_acid):
	prefix = amino_acid[0]
	if prefix == 'p':
		if amino_acid in phosphorylation_list:
			return "phosphorylate "
		else:
			print("The only amino acids that support phosphorylation are S, T and Y")
			sys.exit(0)
	elif prefix == 'n':
		if amino_acid in s_nitrosylation_list: # Keep in mind that we will need an elif here for the rest of n's
			return "nitrosylate "
		elif amino_acid in p_hydroxylation_list:
			return "proline-hydroxylate "
		elif amino_acid in nitration_list:
			return "nitrate "
		else:
			print("The PTMs that have the n as prefix are s-nitrosylation and p-hydroxylation (maybe nitration also). For these PTMs, the only supported amino acids are C, and P (maybe W and Y)")
			sys.exit(0)
	elif prefix == 'c':
		if amino_acid in citrullination_list: # Keep in mind that we will need an elif here for the rest of c's
			return "citrullinate "
		elif amino_acid in carbamylation_list:
			return "carbamylate "
		else:
			print("The PTMs that have the c as prefix are carbamylation and citrullination (maybe c_oxidation also). For these PTMs, the only supported amino acids are C, K and R")
			sys.exit(0)
	elif prefix == 'a':
		if amino_acid in acetylation_list:
			return "acetylate "
		else:
			print("The PTM that has the a as prefix is acetylation. For these PTM, the only supported amino acid is K.")
			sys.exit(0)
	elif prefix == 'm':
		if amino_acid in methylation_list:
			if amino_acid[1] == 'm':
				return "methylate "
			elif amino_acid[1] == 'd':
				return "di-methylate "
			elif amino_acid[1] == 't':
				return "tri-methylate "
			else:
				print("PTM chosen is methylation. You can only mono-methylate (mm), di-methylate (md) or tri-methylate (mt).")
				sys.exit(0)
		else:
			print("PTM chosen is methylation. You can only mono-methylate (mm), di-methylate (md) or tri-methylate (mt).")
			sys.exit(0)
	elif prefix == 'x':
		if amino_acid in c_oxidation_list:
			if amino_acid[1] == 'h':
				return "cysteine-hydroxydate "
			elif amino_acid[1] == 'o':
				return "cysteine-oxydate "
			elif amino_acid[1] == 'd':
				return "cysteine-dioxydate "
			else:
				print("PTM chosen is cysteine oxidation. You can only cysteine-hydroxidate (xh), cysteine-oxidate (xo) or cysteine-dioxidate (xd).")
				sys.exit(0)
		else:
			print("PTM chosen is methylation. You can only mono-methylate (mm), di-methylate (md) or tri-methylate (mt).")
			sys.exit(0)
	elif prefix == 'o':
		if amino_acid in m_oxidation_list:
			return "methionine-oxidization "
		else:
			print("PTM chosen is methionine oxidization. For these PTM, the only supported amino acid is M.")
			sys.exit(0)
	else:
		print("Wrong PTM prefix, check PTM notation")
		sys.exit(0)



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
				 'CYS':[["CA", "C", "CB", "SG", "HG1"]],
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
					   'CYS':[[["CA", "C"], ["CA", "CB"], ["CB", "SG"], ["SG", "HG1"]]],
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



## FEATURE EXTRACTION FOR RF INPUT

def extract_positions(position, peptide, feature_size):
	scale = range(-int(np.floor(feature_size/2)), int(np.floor(feature_size/2)) + 1)
	position_list = []
	append_left = 0
	append_right = 0
	for x in scale:
		if (position + x < 1):
			position_list.append(0)
			append_left += 1
		elif (position + x > len(peptide)):
			position_list.append(0)
			append_right += 1
		else:
			position_list.append(position + x)
	return position_list, append_left, append_right

def extract_features(peptide, MHC, frequencies):

	# Parameters (that more or less give the best RF performance)
	windows_size = 5
	feature_size = 5
	wind = int(np.floor(windows_size/2))

	# Loading the BLOSUM62 matrix
	substitution_matrix = substitution_matrices.load("BLOSUM62")
	substitution_matrix['X'] = substitution_matrix['*'] # X seems to be something else?

	# Extract Ngram subsequences
	aa_neighbors = [''.join(w) for w in list(ngrams(peptide, windows_size, pad_left=True, pad_right=True, 
													left_pad_symbol='X', right_pad_symbol='X'))][wind:]
	aa_neighbors = aa_neighbors[:len(aa_neighbors) - wind]
	
	# Extract all features for each subsequence
	sub_seq_feature_list = []
	for k, sub_seq in enumerate(aa_neighbors):
	
		# Break subsequence down to amino_acids
		pep_sequence = list(sub_seq)
		
		# For each amino-acid, extract frequencies, important positions, and BLOSUM62!
		feature_list = []
		for j, amino_acid in enumerate(pep_sequence):
			
			blosum62_features = np.array(substitution_matrix[amino_acid])
			freq_features = frequencies[(frequencies['allele'] == MHC) & (frequencies['length'] == len(peptide))]
			position_list, append_left, append_right = extract_positions(j + k + 1 - wind, peptide, feature_size)
			freq_features = freq_features[freq_features['position'].isin(position_list)][amino_acid].values
			freq_features = np.hstack([blosum62_features, np.zeros(append_left), freq_features, np.zeros(append_right)])
			feature_list.append(freq_features) 
		
		features = np.array(np.hstack(feature_list))
		posn_feature = (k == 0)
		posnplus1_feature = (k == 1)
		posc_feature = (k == len(peptide) - 1) # Get those out!
		features = np.hstack([features, [posn_feature], [posnplus1_feature], [posc_feature]])
		sub_seq_feature_list.append(features)

	X = np.array(np.array(sub_seq_feature_list, dtype=np.float), dtype=object)

	return X