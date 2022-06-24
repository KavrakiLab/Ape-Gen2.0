import sys
sys.path.insert(1, '/data/pymol_scripts')
import pymol
import pytms

reinitialize
set pdb_retain_ids

input_file = sys.argv[1]
selection = sys.argv[2]
output_file = sys.argv[3]

cmd.load(input_file)
cmd.select(name = "sele", selection = "resi " + selection)
phosphorylate selection="sele", optimize = 1
cmd.save(output_file)

deselect