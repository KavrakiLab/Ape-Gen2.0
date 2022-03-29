import pandas as pd
from biopandas.pdb import PandasPdb
import numpy as np
import sys, importlib

from helper_scripts.Ape_gen_macros import merge_and_tidy_pdb, remove_remarks_and_others_from_pdb

from constraint import *
import itertools

def actual_function():

	# Load flexible receptor
	flex_receptor = PandasPdb()
	flex_receptor.read_pdb("./intermediate_files/1/SMINA_data/flexible_receptors/receptor_58.pdb")
	flex_receptor_df = flex_receptor.df['ATOM']

	# Load rigid receptor
	rigid_receptor = PandasPdb()
	rigid_receptor.read_pdb("./intermediate_files/1/SMINA_data/receptor_for_smina.pdb")
	rigid_receptor_df = rigid_receptor.df['ATOM']
	rigid_receptor_df = rigid_receptor_df[rigid_receptor_df['chain_id'] == 'A'].copy()

	flexible_residues_list = pd.unique(flex_receptor_df['residue_number']).tolist() # Alternative: Load those from object?

	for residue in flexible_residues_list:
		
		flexible_residue = flex_receptor_df[flex_receptor_df['residue_number'] == residue].copy()

		# Create another dataframe with the ATOM names altered:
		rigid_receptor_altered = rigid_receptor_df[rigid_receptor_df['residue_number'] == residue].copy()
		print(rigid_receptor_altered)
		print(flexible_residue)
		input()
		
		# Remove backbone N, H, O
		rigid_receptor_altered = rigid_receptor_altered.iloc[2:, :]
		rigid_receptor_altered = rigid_receptor_altered[rigid_receptor_altered['atom_name'] != 'O']

		# Remove all hydrogens not associated with an N or an O (are these non-polar?)
		altered_receptor_index = 0
		wanted_indexes_list = []
		while altered_receptor_index < rigid_receptor_altered.shape[0]:
			#print(altered_receptor_index)
			altered_atom_name = rigid_receptor_altered.iloc[altered_receptor_index, rigid_receptor_altered.columns.get_loc('element_symbol')]
			if altered_atom_name in ['C', 'S']:
				wanted_indexes_list.append(altered_receptor_index)
			if altered_atom_name in ['N', 'O']:
				wanted_indexes_list.append(altered_receptor_index)
				altered_receptor_index += 1
				while altered_receptor_index < rigid_receptor_altered.shape[0] and rigid_receptor_altered.iloc[altered_receptor_index, rigid_receptor_altered.columns.get_loc('element_symbol')] == 'H':
					wanted_indexes_list.append(altered_receptor_index)
					altered_receptor_index += 1
				altered_receptor_index -= 1
			altered_receptor_index += 1
		
		rigid_receptor_altered = rigid_receptor_altered.iloc[wanted_indexes_list, :].copy()

		atoms_to_iterate = pd.unique(rigid_receptor_altered['element_symbol']).tolist()
		for atom in atoms_to_iterate:
			rigid_receptor_altered_copy = rigid_receptor_altered[rigid_receptor_altered['element_symbol'] == atom].copy()
			flexible_residue_copy = flexible_residue[flexible_residue['element_symbol'] == atom].copy()
			rigid_receptor_df.loc[[x - 1 for x in rigid_receptor_altered_copy['atom_number'].tolist()], 
								  ['x_coord', 'y_coord', 'z_coord']] = flexible_residue_copy[['x_coord', 'y_coord', 'z_coord']].to_numpy()			  

		#altered_receptor_index = 0
		#for index, flexible_row in flexible_residue.iterrows():

		#	new_x_coord = flexible_row['x_coord']
		#	new_y_coord = flexible_row['y_coord']
		#	new_z_coord = flexible_row['z_coord']

		#	altered_entry_found = False
		#	while not altered_entry_found:
		#		altered_receptor_index += 1
		#		altered_atom_name = rigid_receptor_altered.iloc[altered_receptor_index - 1, rigid_receptor_altered.columns.get_loc('element_symbol')]
		#		altered_atom_number = rigid_receptor_altered.iloc[altered_receptor_index - 1, rigid_receptor_altered.columns.get_loc('atom_number')]
		#		if flexible_row['atom_name'] == altered_atom_name:
		#			rigid_receptor_df.loc[rigid_receptor_df['atom_number'] == altered_atom_number, 'x_coord'] = new_x_coord
		#			rigid_receptor_df.loc[rigid_receptor_df['atom_number'] == altered_atom_number, 'y_coord'] = new_y_coord
		#			rigid_receptor_df.loc[rigid_receptor_df['atom_number'] == altered_atom_number, 'z_coord'] = new_z_coord
		#			altered_entry_found = True
		#	print(residue, index, altered_receptor_index)

		print(rigid_receptor_df[rigid_receptor_df['residue_number'] == residue])
		input()

def what_does_the_fox_say():
  print("vixens cry")

def actual_function_2():

	original_ppdb = PandasPdb()
	original_ppdb.read_pdb("./intermediate_files/1/SMINA_data/receptor_for_smina.pdb")
	original_pdb_df = original_ppdb.df['ATOM']
		
	flexible_ppdb = PandasPdb()
	flexible_ppdb.read_pdb("./intermediate_files/1/SMINA_data/flexible_receptors/receptor_47.pdb")
	flexible_pdb_df = flexible_ppdb.df['ATOM']

	flexible_residues = pd.unique(flexible_pdb_df['residue_number'])
	original_ppdb.df['ATOM'] = original_pdb_df[(~(original_pdb_df['residue_number'].isin(flexible_residues))) | (original_pdb_df['atom_name'].isin(["N", "O"]))]
	original_ppdb.to_pdb(path="/data/temp.pdb", records=['ATOM'], gz=False, append_newline=True)
	flexible_ppdb.to_pdb(path="/data/temp2.pdb", records=['ATOM'], gz=False, append_newline=True)
	merge_and_tidy_pdb(["/data/temp.pdb", "/data/temp2.pdb"], "/data/temp3.pdb")

def renaming_function():

	# Making the CONECT list first:
	edge_list = []
	taken = remove_remarks_and_others_from_pdb("./intermediate_files/1/SMINA_data/flexible_receptors/receptor_24.pdb", 
																		           records=("CONECT "))
	conect_fields = ''.join(taken)
	for line in conect_fields.split('\n'):
		cleaned_line = [elem for elem in line.strip().split(' ') if len(elem) > 0]
		try:
			pivot = cleaned_line[1]
			for atom_index in cleaned_line[2:]:
				edge_list.append([int(pivot), int(atom_index)])
		except IndexError:
			pass

	original_ppdb = PandasPdb()
	original_ppdb.read_pdb("./intermediate_files/1/SMINA_data/receptor_for_smina.pdb")
	original_pdb_df = original_ppdb.df['ATOM']

	flexible_ppdb = PandasPdb()
	flexible_ppdb.read_pdb("./intermediate_files/1/SMINA_data/flexible_receptors/receptor_24.pdb")
	flexible_pdb_df = flexible_ppdb.df['ATOM']

	flexible_residues = pd.unique(flexible_pdb_df['residue_number'])
	list_of_dataframes = []
	for flex_residue in flexible_residues:
		sub_origin_pdb = original_pdb_df[(original_pdb_df['residue_number'] == flex_residue) & (original_pdb_df['chain_id'] == 'A')].copy()
		sub_pdb = flexible_pdb_df[flexible_pdb_df['residue_number'] == flex_residue].copy()
		atom_indexes = sub_pdb['atom_number'].tolist()
		sub_edge_list = [elem for elem in edge_list if ((elem[0] in atom_indexes) and (elem[1] in atom_indexes))]
		residue = (sub_pdb['residue_name'].tolist())[0]

		# Remember that CA and C are in the backbone, and are not moving at all, so their co-ordinates will be the same
		# CA location
		CA_coords = np.array(sub_origin_pdb[sub_origin_pdb['atom_name'] == 'CA'][['x_coord', 'y_coord', 'z_coord']])
		loc_indexes = np.all(np.isclose(CA_coords, np.array(sub_pdb[['x_coord', 'y_coord', 'z_coord']]), 
			                              rtol=1e-05, atol=1e-08, equal_nan=False), axis = 1)
		CA_loc = (sub_pdb.loc[loc_indexes, 'atom_number'].values)[0]

		# C location
		C_coords = np.array(sub_origin_pdb[sub_origin_pdb['atom_name'] == 'C'][['x_coord', 'y_coord', 'z_coord']])
		loc_indexes = np.all(np.isclose(C_coords, np.array(sub_pdb[['x_coord', 'y_coord', 'z_coord']]), 
			                              rtol=1e-05, atol=1e-08, equal_nan=False), axis = 1)
		C_loc = (sub_pdb.loc[loc_indexes, 'atom_number'].values)[0]
		
		matching = csp_solver(sub_edge_list, residue, atom_indexes, CA_loc, C_loc)

		sub_pdb = sub_pdb.drop(columns='atom_name').merge(matching, how='inner', on='atom_number')
		list_of_dataframes.append(sub_pdb)
	
	flexible_ppdb.df['ATOM'] = pd.concat(list_of_dataframes)

	original_ppdb.df['ATOM'] = original_pdb_df[(~(original_pdb_df['residue_number'].isin(flexible_residues))) | (original_pdb_df['atom_name'].isin(["N", "O", "H"]))]
	original_ppdb.to_pdb(path=filestore + "/SMINA_data/temp.pdb", records=['ATOM'], gz=False, append_newline=True)
	flexible_ppdb.to_pdb(path=filestore + "/SMINA_data/flexible_receptors/receptor_" + str(peptide_index) + ".pdb", records=['ATOM'], gz=False, append_newline=True)
	merge_and_tidy_pdb([filestore + "/SMINA_data/temp.pdb", 
											filestore + "/SMINA_data/flexible_receptors/receptor_" + str(peptide_index) + ".pdb"],
											minimized_receptor_loc)
	os.remove(filestore + "/SMINA_data/temp.pdb")

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
	print(residue)
	solution = problem.getSolutions()[0]
	return pd.DataFrame(data={'atom_name': solution.keys(), 'atom_number': list(solution.values())})
