import os
import shutil
import sys
import re
from pathlib import Path

import pandas as pd

from pdbtools import pdb_merge, pdb_tidy, pdb_reatom, pdb_sort, pdb_selchain

from constraint import *

## MACROS

# Three-to-one (and vice versa) AA transformation
standard_three_to_one_letter_code = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', \
		                  			 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', \
		                  			 'GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', 'ILE':'I', \
		                  			 'LEU':'L', 'MET':'M', 'PHE':'F', 'TYR':'Y', 'TRP':'W'}
standard_one_to_three_letter_code = {value: key for key, value in standard_three_to_one_letter_code.items()}


non_standard_three_to_one_letter_code = {'SEP':'pS', 'TPO':'pT', 'PTR':'pY'}
non_standard_one_to_three_letter_code = {value: key for key, value in non_standard_three_to_one_letter_code.items()}

all_three_to_one_letter_codes = dict(standard_three_to_one_letter_code, **non_standard_three_to_one_letter_code)
all_one_to_three_letter_codes = dict(standard_one_to_three_letter_code, **non_standard_one_to_three_letter_code)


## FUNCTIONS

def AA_error_checking(amino_acid):
	if amino_acid not in standard_three_to_one_letter_code.values():
		print("The provided amino acid in the sequence is wrong")
		sys.exit(0)

def initialize_dir(dirname):
	if os.path.exists(dirname):
		shutil.rmtree(dirname)
	else:
		pass
	os.makedirs(dirname)

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

def merge_and_tidy_pdb(list_of_pdbs, dst):
	merged = pdb_merge.run(pdb_merge.check_input(list_of_pdbs))
	sorteded = pdb_sort.run(merged, sorting_keys='-RC') # Potentially breaking -> Not working?
	tidied = pdb_tidy.run(sorteded, strict=True)
	reatomed = pdb_reatom.run(tidied, starting_value=1)
	with open(dst, 'w') as pdb_file:
		pdb_file.write(''.join(reatomed))
	pdb_file.close()

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

def delete_elements(pdb_file, element_set):
	records = ('ATOM', 'HETATM', 'ANISOU')
	fhandle = open(pdb_file, 'r')
	for line in fhandle:
		if line.startswith(records):
			if line[13:16].strip() in element_set:
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
acetylation_list = ['aK'] # Check details on pytms
carbamylation_list = ['cK'] # Check details on pytms
citrullination_list = ['cR'] 
methylation_list = ['m1K', 'm2K', 'm3K'] # Check details on pytms
nitration_list = ['nY', 'nW'] # Check details on pytms
s_nitrosylation_list = ['nC']
p_hydroxylation_list = ['n4rP', 'n4sP'] # Check details on pytms
malondialdehyde_list = ['maK'] # Check details on pytms
c_oxidation_list = ['csoC', 'csxC'] # Check details on pytms
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
		else:
			print("The PTMs that have the n as prefix are s-nitrosylation and nitration_list (maybe p-hydroxylation also). For these PTMs, the only supported amino acids are C, P, W and Y")
			sys.exit(0)
	elif prefix == 'c':
		if amino_acid in citrullination_list: # Keep in mind that we will need an elif here for the rest of c's
			return "citrullinate "
		else:
			print("The PTMs that have the c as prefix are carbamylation and citrullination (maybe c_oxidation also). For these PTMs, the only supported amino acids are C, K and R")
			sys.exit(0)
	else:
		print("Wrong PTM prefix, check PTM notation")
		sys.exit(0)



## RESIDUE RENAMING AFTER SMINA FLEXIBILITY OUTPUT

def csp_solver(edge_list, residue, atom_indexes, CA_loc, C_loc):

	# Note to change the 4-letter atoms if need be!
	atom_dict = {'ALA':["CA", "C", "CB"],
				 'VAL':["CA", "C", "CB", "CG1", "CG2"],
				 'ILE':["CA", "C", "CB", "CG1", "CG2", "CD"],
				 'LEU':["CA", "C", "CB", "CG", "CD1", "CD2"],
				 'MET':["CA", "C", "CB", "CG", "SD", "CE"],
			     'PHE':["CA", "C", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
				 'TYR':["CA", "C", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "HH"],
				 'TRP':["CA", "C", "CB", "CG", "CD1", "CD2", "NE1", "HE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
				 'SER':["CA", "C", "CB", "OG", "HG1"],
				 'THR':["CA", "C", "CB", "OG1", "HG1", "CG2"],
				 'ASN':["CA", "C", "CB", "CG", "OD1", "ND2", "HD21", "HD22"],
				 'GLN':["CA", "C", "CB", "CG", "CD", "OE1", "NE2", "HE21", "HE22"],
				 'CYS':["CA", "C", "CB", "SG"],
				 'GLY':["CA", "C"],
				 'PRO':["CA", "C", "CB", "CG", "CD"],
				 'ARG':["CA", "C", "CB", "CG", "CD", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"],
				 'HIS':["CA", "C", "CB", "CG", "ND1", "HD1", "CE1", "NE2", "CD2"],
				 'LYS':["CA", "C", "CB", "CG", "CD", "CE", "NZ", "HZ1", "HZ2", "HZ3"],
				 'ASP':["CA", "C", "CB", "CG", "OD1", "OD2"],
				 'GLU':["CA", "C", "CB", "CG", "CD", "OE1", "OE2"]
			    }

	constraint_dict = {'ALA':[["CA", "C"], ["CA", "CB"]],
					   'VAL':[["CA", "C"], ["CA", "CB"], ["CB", "CG1"], ["CB", "CG2"]],
					   'ILE':[["CA", "C"], ["CA", "CB"], ["CB", "CG1"], ["CB", "CG2"], ["CG1", "CD"]],
					   'LEU':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD1"], ["CG", "CD2"]],
					   'MET':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "SD"], ["SD", "CE"]],
					   'PHE':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD1"], ["CG", "CD2"], ["CD1", "CE1"], ["CD2", "CE2"], ["CE1", "CZ"], ["CE2", "CZ"]],
					   'TYR':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD1"], ["CG", "CD2"], ["CD1", "CE1"], ["CD2", "CE2"], ["CE1", "CZ"], ["CE2", "CZ"], ["CZ", "OH"], ["OH", "HH"]],
					   'TRP':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD1"], ["CG", "CD2"], ["CD1", "NE1"], ["NE1", "CE2"], ["NE1", "HE1"], ["CD2", "CE2"], ["CD2", "CE3"], ["CE2", "CZ2"], ["CE3", "CZ3"], ["CZ2", "CH2"], ["CH2", "CZ3"]],
					   'SER':[["CA", "C"], ["CA", "CB"], ["CB", "OG"], ["OG", "HG1"]],
					   'THR':[["CA", "C"], ["CA", "CB"], ["CB", "OG1"], ["OG1", "HG1"], ["CB", "CG2"]],
					   'ASN':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "ND2"], ["ND2", "HD21"], ["ND2", "HD22"]],
					   'GLN':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "OE1"], ["CD", "NE2"], ["NE2", "HE21"], ["NE2", "HE22"]],
					   'CYS':[["CA", "C"], ["CA", "CB"], ["CB", "SG"], ["SG", "HG1"]],
					   'GLY':[["CA", "C"]],
					   'PRO':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"]],
					   'ARG':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "NE"], ["NE", "HE"], ["NE", "CZ"], ["CZ", "NH1"], ["CZ", "NH2"], ["NH1", "HH11"], ["NH1", "HH12"], ["NH2", "HH21"], ["NH2", "HH22"]],
					   'HIS':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "ND1"], ["ND1", "CE1"], ["ND1", "HD1"], ["CE1", "NE2"], ["NE2", "CD2"], ["CD2", "CG"]],
					   'LYS':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "CE"], ["CE", "NZ"], ["NZ", "HZ1"], ["NZ", "HZ2"], ["NZ", "HZ3"]],
					   'ASP':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "OD1"], ["CG", "OD2"]],
					   'GLU':[["CA", "C"], ["CA", "CB"], ["CB", "CG"], ["CG", "CD"], ["CD", "OE1"], ["CD", "OE2"]]
					   }

	# Formulating the atom matching problem as a CSP:
	problem = Problem()

	# Adding atoms as variables to the CSP (along with possible index values they can take (CA and C are known):
	problem.addVariable("CA", [CA_loc])
	problem.addVariable("C", [C_loc])
	problem.addVariables((atom_dict[residue])[2:], atom_indexes)

	# Adding Constraints to the CSP
	# 1. No atoms share the same index
	problem.addConstraint(AllDifferentConstraint())

	# 2. Topological constraints that are associated with each residue:
	for elem in constraint_dict[residue]:
		problem.addConstraint(lambda a, b: ([a,b] in edge_list) and ([b,a] in edge_list), elem)

	# 3. Find problem solution
	solution = problem.getSolutions()[0]
	return pd.DataFrame(data={'atom_name': solution.keys(), 'atom_number': list(solution.values())})

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