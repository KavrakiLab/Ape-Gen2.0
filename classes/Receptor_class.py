from helper_scripts.Ape_gen_macros import merge_and_tidy_pdb

from biopandas.pdb import PandasPdb
import pandas as pd

from subprocess import call
import sys
import os

# PDBFIXER
from pdbfixer import PDBFixer
from openmm.app import PDBFile

class Receptor(object):

	def __init__(self, allotype, pdb_filename):
		self.pdb_filename = pdb_filename
		self.pdbqt_filename = None
		self.allotype = allotype
		self.flexible_residues = self.load_flexible_residues()
		self.doMinimization = True
		self.useSMINA = True

	@classmethod
	def frompdb(cls, pdb_filename): # Checking if .pdb file is ok maybe?
		return cls(allotype = "In PDB", pdb_filename = receptor_class)

	@classmethod
	def fromredock(cls, peptide_input):
		return cls(allotype = "REDOCK", pdb_filename = peptide_input)

	@classmethod
	def fromallotype(cls, allotype):

		# Check #1: Existing structures
		templates = pd.read_csv("./template_files/receptor-class-templates.csv")
		if(allotype in templates['Allotype'].tolist()):
			print("Allotype found in our structural DB!")
			pdb_filename = templates[templates['Allotype'] == allotype]['Template'].values[0]
			return(cls(allotype = allotype, pdb_filename = './templates/' + pdb_filename))
		
		# Check #2: Supertype representatives
		print("Allotype not found in our structural DB. Will resort to modelling the structure...")
		print(" Homology modelling (MODELLER) is used. Trying to fetch supertype representative as the template:")
		supertypes = pd.read_csv("./template_files/supertype_templates.csv")
		if(allotype in supertypes['Allotype'].tolist()):
			print("Supertype for given allotype found in our DB!")
			receptor_template = templates[templates['Allotype'] == allotype]['Template'].values[0]
			print("Using " + receptor_class + "as a representative template for " + allotype)
			
		else:
			print("Warning: Allele cannot be found in the supertype database ... Defaulting to 2v2w.pdb")
			receptor_template = "./templates/2v2w/pdb"
		##MODELLER FUNCTION CALL
		return(cls(allotype = allotype, pdb_filename = pdb_filename))

	@classmethod
	def fromfasta(cls, sequence):
		#Idea is to call the sequence similarity to fetch allele -> fetch representative -> MODELLER	
		return cls(allotype = None, pdb_filename = None)

	@staticmethod	
	def load_flexible_residues():
		file = open('./template_files/flex_res.txt', "r")
		flexible_residues = file.readline().strip()
		return flexible_residues

	def add_sidechains(self, filestore):
		fixer = PDBFixer(filename=self.pdb_filename)
		fixer.findMissingResidues()
		fixer.removeHeterogens(True) #  True keeps water molecules while removing all other heterogens, REVISIT!
		fixer.findMissingAtoms()
		fixer.addMissingAtoms()
		fixer.addMissingHydrogens(7.0) # Ask Mauricio about those
		#fixer.addSolvent(fixer.topology.getUnitCellDimensions()) # Ask Mauricio about those
		PDBFile.writeFile(fixer.topology, fixer.positions, open(self.pdb_filename, 'w'))

	def prepare_for_scoring(self, filestore):

		prep_receptor_loc = "/conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
		pdbqt_to_pdb_loc = "/conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py"
		
		self.pdbqt_filename = filestore + "/SMINA_data/receptor_for_smina.pdbqt"
		call(["python2.7 " + prep_receptor_loc + " -r " + self.pdb_filename + " -o " + self.pdbqt_filename + " -A None -U lps > " + filestore + "/SMINA_data/prepare_receptor4.log 2>&1"], shell=True)
		call(["python2.7 " + pdbqt_to_pdb_loc + " -f " + self.pdbqt_filename + " -o " + filestore + "/SMINA_data/receptor_for_smina_temp.pdb > " + filestore + "/SMINA_data/pdbqt_to_pdb.log 2>&1"], shell=True)

		# Before we continue here, an issue seems to arise. pdbqt_to_pdb.py introduces some segment identifiers that need to be removed?
		self.pdb_filename = filestore + "/SMINA_data/receptor_for_smina_temp.pdb"
		ppdb_receptor = PandasPdb()
		ppdb_receptor.read_pdb(self.pdb_filename)
		pdb_df_receptor = ppdb_receptor.df['ATOM']
		ppdb_receptor.df['ATOM']['segment_id'] = ''
		ppdb_receptor.to_pdb(path=self.pdb_filename, records=None, gz=False, append_newline=True)

		# Adding the following lines to properly have TER and END fields (hence the temp file here, maybe there's a better way to do this)
		self.pdb_filename = filestore + "/SMINA_data/receptor_for_smina.pdb"
		merge_and_tidy_pdb([filestore + "/SMINA_data/receptor_for_smina_temp.pdb"], self.pdb_filename)
		os.remove(filestore + "/SMINA_data/receptor_for_smina_temp.pdb")