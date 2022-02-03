# New Ape-Gen version:

## Installation

For now, only going to `\Docker_Image` and running:
```
docker build . -t kavrakilab/apegen2.0
```

works. After this, just running the .sh file gets you into the container, and Ape-Gen can be run from there. 

An additional step is required here: One needs to go to:

```
cd /conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/
```

and modify the following lines of `prepare_ligand4.py`:

```
if o in ('-l', '--l'):
	#ligand_filename = a
	ligand_filename = os.path.basename(a) 
```

to

```
if o in ('-l', '--l'):
	ligand_filename = a
	#ligand_filename = os.path.basename(a)
```
Not sure why this is happening in autodocktools, will see if I can manually overwrite the file while building the image. 


For now, peptide sequence (along with phosphorylated positions) + HLA allotype works, but native + other types of inputs will follow:

```
python New_APE-Gen.py ARpSEDEVILpS HLA-A*11:01 --debug
```

## TO-DOs:

### Main Workflow:

#### Minor issues:
- Manually overwrite `prepare_ligand4.py` with the changes that do not cause errors
- Get debug message inside the print statements happening in the classes
- Fix `initialize_dir` function to accept a list of directories and avoid calling the function multiple times
- Consider peptide methods which receptor is an input to be transferred to the `pMHC_class`, and a pMHC object defined
- `peptide_index` in methods is annoying, maybe replace it with a new field
- `simtk.openmm` import warning coming from PDBFixer, wait for a new version or modify the .xml: https://github.com/openmm/pdbfixer/issues/233
- See if you can fix the script that replaces the flexible residues to the receptor (`ProDy` installation could be avoided)
- Better folder names for intermediate results (it would be nice if they are numbered so that they are ordered nicely)
- Function documentation needs to be done thouroughly at some point

#### Major issues:
- Implement many rounds (currently not working properly)
- Implement other inputs (native, HLA sequence, etc.)
- Testing/Testing/Testing...

### PTMs:

- Implement other PTMs:
	- Already done:
		- Phosphorylation
		- Citrullination
		- S-Nitrosylation
	- To Do:
		- Acetylation
		- Carbamylation
		- Methylation
		- Nitration
		- P-hydroxylation
		- Malondialdehyde adducts
		- C-Oxidation
		- M-Oxidation
- Ask Mauricio/Dinler about N-termini PTMs
- Bring also Deamidation issue

### OpenMM:

- Re-writing OpenMM and improve the code
- Implement PTMs in the OpenMM step (see phosphorylation example that is working)
- Update forcefield parameters for PTMs with the one released this year from Mackerel's group
	- https://www.charmm.org/charmm/resources/charmm-force-fields/
	- https://github.com/openmm/openmmforcefields/tree/master/charmm#Converting

### Actual Experiments:

- Fetch as many templates as you can and assess experiments below in terms of RMSD?
- Test RCD parameters (mainly number of conformations that are needed)
- Check Ape-gen cycles
- Phosphorylated peptides: What experiments can be done here? Using Anja's scoring function maybe?

### Other possible Expansions:

- GUI
- Ensemble of receptor conformations in RCD step (maybe after that too?)
- AlphaFold / RosettaTFold (wait for an API like thing) (could I also generate a bunch of HLAs offline and assess instead?)
- Possibly insert GNINA in workflow somehow? (and test it's RMSD to vinardo/vina)
- Think about alternative scoring function? )and test it's RMSD to vindardo/vina/CNN)

