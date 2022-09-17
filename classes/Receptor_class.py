from helper_scripts.Ape_gen_macros import remove_file, merge_and_tidy_pdb, copy_file, initialize_dir, verbose

from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np

from subprocess import call
import sys
import time
import re
import os

from Bio import Align
from Bio import SeqIO
from Bio import pairwise2

# MODELLER
try:
	from modeller import *
	from modeller.automodel import *
	modeller_import = None
except:
	modeller_import = ImportError

def model_single_opt(filestore, num_models=10):
	# if the modeller import was unsuccessful, quit
	if (modeller_import == ImportError):
		print("Error with importing Modeller: Make sure license key is correct.")
		sys.exit(0)

	log.none()
	env = Environ()

	# Give less weight to all soft-sphere restraints:
	env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)

	#Considering heteroatoms and waters molecules
	env.io.hetatm = env.io.water = True

	# Directories with input atom files:
	# env.io.atom_files_directory = './:../atom_files'

	# Modelling 'sequence' with file.ali
	a = AutoModel(env, alnfile='target_sequence-receptor_template.ali', knowns='receptor_templateA',
				  sequence='target_sequence', assess_methods=(assess.DOPE, assess.GA341))

	# Generating 3 models
	a.starting_model = 1
	a.ending_model = int(num_models)

	# Very thorough Variable Target Function Method (VTFM) optimization:
	a.library_schedule = autosched.slow
	a.max_var_iterations = 300

	# Thorough MD optimization:
	a.md_level = refine.slow
	 
	# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
	a.repeat_optimization = 2
	a.max_molpdf = 1e6

	a.make()

	# Get a list of all successfully built models from a.outputs
	ok_models = [x for x in a.outputs if x['failure'] is None]

	# Rank the models by DOPE score
	key = 'DOPE score'
	ok_models.sort(key=lambda a: a[key])

	# Get top model
	m = ok_models[0]
	if verbose(): print("Top model: %s (DOPE score %.100f)" % (m['name'], m[key]))

	return m['name'], m[key]

def align_2d(filestore):

	log.none()

	env = Environ()
	aln = Alignment(env)
	mdl = Model(env, file='receptor_template', model_segment=('FIRST:A','LAST:A'))
	aln.append_model(mdl, align_codes='receptor_templateA', atom_files='receptor_template.pdb')
	aln.append(file='target_sequence.pir', align_codes='target_sequence')
	aln.align2d()
	aln.write(file='target_sequence-receptor_template.ali', alignment_format='PIR')
	aln.write(file='target_sequence-receptor_template.pap', alignment_format='PAP')

def model_receptor(allele_sequence, peptide_sequence, filestore):
	# if the modeller import was unsuccessful, quit
	if (modeller_import == ImportError):
		print("Error with importing Modeller: Make sure license key is correct.")
		sys.exit(0)

	# Routine for finding the best sequence match for homology modelling, as well as modelling itself

	# First we search for the closest sequence match
	if verbose(): print("Searching for the closest match in terms of sequence:")
	best_record_list = []
	best_score = 0 
	for seq_record in SeqIO.parse("./helper_files/template_sequences.fasta", "fasta"):

		# Do a quick alignment
		
		alignments = pairwise2.align.globalxx(allele_sequence, str(seq_record.seq))
		if(alignments[0].score == best_score):
			best_record_list.append(seq_record.id)
		if(alignments[0].score > best_score):
			best_score = alignments[0].score
			best_record_list = []
			best_record_list.append(seq_record.id)
	if verbose(): print(best_record_list)
	if verbose(): print("Closest match are the following alleles: ", best_record_list)

	# Secondly, emphasizing more on the peptide binding pocket, we fetch the pdb code that has the peptide
	# with the closest sequence similarity to the peptide to be modelled. This is coming from the hypothesis that
	# this structure will have similar binding cleft to accomodate the peptide?

	# This is repeated code btw, see if you can make a function for this, it would be amazing
	if verbose(): print("Fetching now the most appropriate template that will host the peptide in question:")
	templates = pd.read_csv("./helper_files/Template_Information_notation.csv")
	templates = templates[templates['MHC'].isin(best_record_list)]
	aligner = Align.PairwiseAligner()
	aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
	score_list = []
	template_sequences = templates['peptide'].tolist()
	for template_sequence in template_sequences:
		score_list.append(aligner.score(peptide_sequence, template_sequence))
	templates['peptide_score'] = score_list
	templates = templates[templates['peptide_score'] == templates['peptide_score'].max()].dropna()
	result = templates.sample(n=1)
	pdb_filename = result['pdb_code'].values[0]
	new_allotype = result['MHC'].values[0]
	if verbose(): print("Got " + pdb_filename + "! This template has MHC " + new_allotype)
	copy_file('./new_templates/' + pdb_filename, filestore + '/receptor_template.pdb')

	# Now that we found the pdb that is serving as the HLA template, take the sequence of that template
	# CAREFUL: the sequence is only the A-chain!
	# Solution: Do a global alignment, where you penalize opening gaps, as such, the gap in the beginning is
	# that weird starting aa seq, and ther gap in the end will be the B-chain.
	# THINK ABOUT WHEN THIS FAILS!
	for seq_record in SeqIO.parse("./helper_files/a_chain_sequences.fasta", "fasta"):
		if seq_record.id == pdb_filename:
			a_chain_alignment = pairwise2.align.globalxs(allele_sequence, str(seq_record.seq), open=-.5, extend=0)
	original_sequence = list(str(a_chain_alignment[0][0]))
	gapped_sequence = list(str(a_chain_alignment[0][1]))
	filtered_sequence = [original_sequence[i] for i, x in enumerate(gapped_sequence) if x != '-']
	filtered_sequence = "".join(filtered_sequence) + '*'
	
	wd_to_return_to = os.getcwd()
	os.chdir(wd_to_return_to + '/' + filestore)

	if verbose(): print("Preparing target sequence")
	f = open("target_sequence.pir", 'w')
	f.write(">P1;target_sequence\n")
	f.write("sequence::	 : :	 : :::-1.00:-1.00\n")
	f.write(filtered_sequence)
	f.close()

	align_2d(filestore)

	if verbose(): print("Creating model")
	starttime = time.time()	
	best_model, best_score = model_single_opt(filestore, num_models=2)
	endtime = time.time()
	if verbose(): print("Homology modelling took " + str(endtime - starttime) + " seconds.")
	os.chdir(wd_to_return_to)

	return filestore + '/' + best_model, new_allotype

class Receptor(object):

	def __init__(self, allotype, pdb_filename):
		self.pdb_filename = pdb_filename
		self.pdbqt_filename = None
		self.allotype = allotype
		self.flexible_residues = self.load_flexible_residues()
		self.doMinimization = True
		self.useSMINA = True

	def init_receptor(receptor_class, file_storage, peptide_input, cv=''):
		if verbose(): print("\nProcessing Receptor Input: " + receptor_class)

		if receptor_class.endswith(".pdb"):
			# If the file is .pdb, this will be your template! ##MUST CHECK VALIDITY IN THE FUNCTION
			receptor = Receptor.frompdb(receptor_class)
			receptor_template_file = receptor_class
		elif receptor_class.endswith(".fasta"):
			# If this is a sequence, the template is taken by MODELLER
			initialize_dir(file_storage + '/MODELLER_output')
			receptor = Receptor.fromfasta(receptor_class, peptide_input, file_storage)
			receptor_template_file = receptor.pdb_filename
		elif receptor_class == "REDOCK":
			# If REDOCK, the receptor template is the peptide template!
			receptor = Receptor.fromredock(peptide_input)
			receptor_template_file = peptide.pdb_filename
		else:
			# If this is an allotype specification, fetch template like the peptide!
			initialize_dir(file_storage + '/MODELLER_output')
			receptor = Receptor.fromallotype(receptor_class, peptide_input, file_storage, cv)
			receptor_template_file = receptor.pdb_filename
		return receptor, receptor_template_file

	@classmethod
	def frompdb(cls, pdb_filename): # Checking if .pdb file is ok maybe?
		return cls(allotype="In PDB", pdb_filename=pdb_filename)

	@classmethod
	def fromredock(cls, peptide_input):
		return cls(allotype="REDOCK", pdb_filename=peptide_input)

	@classmethod
	def fromallotype(cls, allotype, peptide_sequence, filestore, cv=''):

		# Check #1: Existing structures
		templates = pd.read_csv("./helper_files/Updated_template_information.csv")

		if cv != '': templates = templates[~templates['pdb_code'].str.contains(cv, case=False)]

		if(allotype in templates['MHC'].tolist()):

			if verbose(): print("Allotype found in our structural DB!")
			templates = templates[templates['MHC'] == allotype]

			if verbose(): print("Will try to get the one that is closer to the peptide_input:")
			# select the one closer to the whole sequence
			aligner = Align.PairwiseAligner()
			aligner.open_gap_score = -0.5
			aligner.extend_gap_score = -0.1
			aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
			score_list = []
			template_sequences = templates['peptide'].tolist()
			for template_sequence in template_sequences:
				score_list.append(aligner.score(peptide_sequence, template_sequence))
			templates['peptide_score'] = score_list
			templates = templates[templates['peptide_score'] == templates['peptide_score'].max()].dropna()
			templates = templates.sample(n=1)
			pdb_filename = templates['pdb_code'].values[0]
			if verbose(): 
				print("\tGot " + pdb_filename + "!")
				print("\tPeptide is " + templates['peptide'].values[0])
				print("\tMHC is " + templates['MHC'].values[0])
			return(cls(allotype=allotype, pdb_filename='./new_templates/' + pdb_filename))
		
		# Check #2: Existing sequence
		if verbose(): print("Allotype not found in our structural DB. Let's see if it's in our sequence DB...")
		for seq_record in SeqIO.parse("./helper_files/MHC_data.fasta", "fasta"):
			if seq_record.id == allotype:
				if verbose(): print("Allotype found in our sequence DB! Modelling it through homology modelling:")
				pdb_filename, new_allotype = model_receptor(str(seq_record.seq), peptide_sequence, filestore)
				return(cls(allotype=new_allotype, pdb_filename=pdb_filename))

		print("Allotype not found in our sequence DB... Please check the list of supported allotypes (or pass down a fasta sequence instead.")
		print("Aborting....")
		sys.exit(1)

	@classmethod
	def fromfasta(cls, fasta_file, peptide_sequence, filestore):

		# If multiple sequences appear in the fasta file, only the last one will be taken into account!
		# The .fasta file location must also be in relation to the main directory. 

		for seq_record in SeqIO.parse(fasta_file, "fasta"):
			sequence = seq_record.seq
		pdb_filename, new_allotype = model_receptor(str(sequence), peptide_sequence, filestore)
		return(cls(allotype=new_allotype, pdb_filename=pdb_filename))

	@staticmethod	
	def load_flexible_residues():
		file = open('./template_files/flex_res.txt', "r")
		flexible_residues = file.readline().strip()
		return flexible_residues

	def prepare_for_scoring(self, filestore, addH, index=""):

		prep_receptor_loc = "/conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
		pdbqt_to_pdb_loc = "/conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py"
		self.pdbqt_filename = filestore + "/receptor_for_smina" + index + ".pdbqt"

		clean = "lps" if addH == "all" else "nphs_lps"
		call(["python2.7 " + prep_receptor_loc + " -r " + self.pdb_filename + " -o " + self.pdbqt_filename + " -A None -U" + clean + " > " + filestore + "/prepare_receptor4.log 2>&1"], shell=True)
		call(["python2.7 " + pdbqt_to_pdb_loc + " -f " + self.pdbqt_filename + " -o " + filestore + "/receptor_for_smina_temp" + index + ".pdb > " + filestore + "/pdbqt_to_pdb.log 2>&1"], shell=True)

		# Before we continue here, an issue seems to arise. pdbqt_to_pdb.py introduces some segment identifiers that need to be removed?
		self.pdb_filename = filestore + "/receptor_for_smina_temp" + index + ".pdb"
		ppdb_receptor = PandasPdb()
		ppdb_receptor.read_pdb(self.pdb_filename)
		pdb_df_receptor = ppdb_receptor.df['ATOM']
		ppdb_receptor.df['ATOM']['segment_id'] = ''
		ppdb_receptor.to_pdb(path=self.pdb_filename, records=None, gz=False, append_newline=True)

		# Adding the following lines to properly have TER and END fields (hence the temp file here, maybe there's a better way to do this)
		self.pdb_filename = filestore + "/receptor_for_smina" + index + ".pdb"
		merge_and_tidy_pdb([filestore + "/receptor_for_smina_temp" + index + ".pdb"], self.pdb_filename)
		remove_file(filestore + "/receptor_for_smina_temp" + index + ".pdb")