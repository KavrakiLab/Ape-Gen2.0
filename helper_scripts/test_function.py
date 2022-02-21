import pandas as pd
from biopandas.pdb import PandasPdb
import numpy as np
import sys, importlib

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