import os
import shutil
import sys
import re
from pathlib import Path

import pandas as pd

from pdbtools import pdb_merge, pdb_tidy, pdb_reatom

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
	tidied = pdb_tidy.run(merged, strict=True)
	reatomed = pdb_reatom.run(tidied, starting_value=1)
	with open(dst, 'w') as anchored_MHC_file:
		anchored_MHC_file.write(''.join(reatomed))
	anchored_MHC_file.close()

def remove_remarks_and_others_from_pdb(pdb_file):
	records = ('ATOM', 'TER', 'END ')
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