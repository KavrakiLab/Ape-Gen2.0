import argparse

def APE_Gen_parser():

	parser = argparse.ArgumentParser(description="Anchored Peptide-MHC Ensemble Generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('peptide_input', type=str, nargs=1, help='Sequence of peptide to dock or pdbfile of crystal structure')
	parser.add_argument('receptor_class', type=str, nargs=1, help='Class descriptor of MHC receptor. Use REDOCK along with crystal input to perform redocking. Or pass a PDB file with receptor')
	parser.add_argument("-n", "--num_cores", type=int, default=8, help='Number of cores to use for RCD and smina computations.')
	parser.add_argument("--num_generated_loops", type=int, default=5000, help='Number of loops to generate with RCD')
	parser.add_argument("--num_loops_for_optimization", type=int, default=100, help='Number of loops to optimize (that will pass as a result of a loop scoring function)')
	parser.add_argument("-t", "--RCD_dist_tol", type=float, default=1.0, help='RCD tolerance (in angstroms) of inner residues when performing IK')
	parser.add_argument("-r", "--rigid_receptor", action="store_true", help='Disable sampling of receptor degrees of freedom specified in flex_res.txt')
	parser.add_argument("-v", "--verbose", action="store_true", help='Print extra information for debugging')
	parser.add_argument("-p", "--save_only_pep_confs", action="store_true", help='Disable saving full conformations (peptide and MHC)')
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
	parser.add_argument("--anchor_selection", type=str, default='secondary', choices=['primary', 'secondary', 'none'], help="Give what type of anchors should be considered in the anchor tolerance step (choose 'primary', 'secondary' or 'none' to skip the anchor tolerance step altogether)")
	parser.add_argument("--cv", type=str, default='', help='ONLY FOR TESTING (to be removed in the final version)')
	parser.add_argument("--loop_score", type=str, default='ICOSA', choices=['RMSD', 'KORP', 'ICOSA', 'none'], help='Choose scoring function for RCD loop scoring (none will avoid scoring altogether)')
	parser.add_argument("--sampling_ratio", type=float, default=0.8, help='The percentage of overall peptide conformations processed (defined by --num_loops_for_optimization flag) that will be coming from RCD sampling.')
	return parser