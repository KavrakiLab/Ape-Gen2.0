import argparse

def APE_Gen_parser():

	parser = argparse.ArgumentParser(description="Anchored Peptide-MHC Ensemble Generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('peptide_input', type=str, nargs=1, help='Sequence of peptide to dock or pdbfile of crystal structure')
	parser.add_argument('receptor_class', type=str, nargs=1, help='Class descriptor of MHC receptor. Use REDOCK along with crystal input to perform redocking. Or pass a PDB file with receptor')
	parser.add_argument("-n", "--num_cores", type=int, default=8, help='Number of cores to use for RCD and smina computations.')
	parser.add_argument("-l", "--num_loops", type=int, default=100, help='Number of loops to generate with RCD. (Note that the final number of sampled conformations may be less due to steric clashes.')
	parser.add_argument("-t", "--RCD_dist_tol", type=float, default=1.0, help='RCD tolerance (in angstroms) of inner residues when performing IK')
	parser.add_argument("-r", "--rigid_receptor", action="store_true", help='Disable sampling of receptor degrees of freedom specified in flex_res.txt')
	parser.add_argument("-v", "--verbose", action="store_true", help='Print extra information for debugging')
	parser.add_argument("-p", "--save_only_pep_confs", action="store_true", help='Disable saving full conformations (peptide and MHC)')
	#parser.add_argument("--n-mer-templates", default="", help='File with n-mer pdb templates.')
	#parser.add_argument("--receptor-class-templates", default="", help='File with pdb receptor class templates')
	#parser.add_argument("--flex_res", default="", help='File with flexible residues')
	parser.add_argument('--anchors', default="", help='delimited list input', type=str)
	parser.add_argument("-a", "--anchor_tol", type=float, default=2.0, help='Anchor tolerance (in angstroms) of first and last backbone atoms of peptide when filtering')
	parser.add_argument("-o", "--score_with_openmm", action="store_true", help='Rescore full conformations with openmm (AMBER)')
	parser.add_argument("-g", "--num_rounds", type=int, default=1, help='Number of rounds to perform.')
	parser.add_argument("-b", "--pass_type", type=str, default='receptor_only', choices=['receptor_only', 'pep_and_recept'], help="When using multiple rounds, pass best scoring conformation across different rounds (choose either 'receptor_only' or 'pep_and_recept')")
	parser.add_argument("-s", "--min_with_smina", action="store_true", help='Minimize with SMINA instead of the default Vinardo')
	parser.add_argument("--use_gpu", action="store_true", help='Use GPU for OpenMM Minimization step')
	parser.add_argument("--clean_rcd", action="store_true", help='Remove RCD folder at the end of each round')
	parser.add_argument("--dir", type=str, default='intermediate_files', help='Location for all the intermediate files')
	parser.add_argument("--force_restart", action="store_true", help='Force restart of APE-Gen *ONLY* in the first round and *ONLY* if no conformations are produced')
	return parser