# New Ape-Gen version:

## Installation

For now, only going to `\Docker_Image` and running:
```
docker build . -t kavrakilab/apegen2.0
```

works. After this, just running the .sh file gets you into the container, and Ape-Gen can be run from there. 

Two additional step are required here: One needs to go to:

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

Secondly, to give a valid MODELLER key, one needs to go to:

```
cd /conda/envs/apegen/lib/modeller-10.2/modlib/modeller/
```

and modify the lines of `config.py` with a valid MODELLER key. Contact the lab to get a valid key.

Finally, download the files from this link:

```
https://rice.box.com/s/q4bh18eqvmvpzulqyfoszebo9jfkivc3
```

and put them in the helper files folder, as these are mandatory for the anchor prediction step!


For now, peptide sequence (along with phosphorylated positions) + HLA allotype works, but native + other types of inputs will follow:

```
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --debug
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --debug --score_with_openmm
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:01 --debug --score_with_openmm
```

## TO-DOs:

### Main Workflow:

#### Minor issues:
- Implement different hydrogen configurations for APE-GEN, like a) Use all hydrogens (new version), b) Only polar hydrogens, c) No hydrogens (old version)
- Maybe have an argument for `autobox` size for SMINA scoring (larger values will result in greater time but much more accurate docking)
- Possibly insert GNINA in workflow somehow? (and test it's RMSD to vinardo/vina)
- Locate all the PTM removal regex expressions and make a function instead
- Make homology modelling into a different script in order not to ask for MODELLER keys for cases where HLA allotype has a known structure.
- Double check on how the flexible residue code assignment works, because it may attempt to enter the `minimized_receptors` folder without a successfull CSP assignment.
- Make `anchor tol` step optional.
- Add `num_of_rounds` in homology modelling as an argument. 
- Make function in the places where there is repeated code (there is one on peptide similarity tempalte selection)
- Have a `write_string` function for all the yields/joins/opens/writes in the code (or maybe for each function do this inside it, the file is overwritten successfully without defining a new one)
- Add `number_of_tries` parameter to OpenMM (maybe also add the Langevin integrator patameters?)
- CSP routine that builds the correct atom names in the flexible file could probably be more optimized (CSP could potentially solve the whole thing instead of per residue?)
- Implement `PDBFixer` steps in a better way (it's the same function in all of the 3 classes basically)
- Make `os.remove` calls a macro function
- The `receptor_template_file` variable in Section 1b is questionable, investigate more
- Manually overwrite `prepare_ligand4.py` with the changes that do not cause errors
- Get debug message inside the print statements happening in the classes
- Fix `initialize_dir` function to accept a list of directories and avoid calling the function multiple times
- Consider peptide methods which receptor is an input to be transferred to the `pMHC_class`, and a pMHC object defined
- `peptide_index` in methods is annoying, maybe replace it with a new class field
- `simtk.openmm` import warning coming from PDBFixer, wait for a new version or modify the .xml: https://github.com/openmm/pdbfixer/issues/233
- Better folder names for intermediate results (it would be nice if they are numbered so that they are ordered nicely)
- Function documentation needs to be done thouroughly at some point
- Thorough input checking (example is peptide sequence in HLA peptide fetching must be an amino acid sequence)

#### Major issues:
- `5TRZ.pdb` is the only peptide with non-canonical anchors in both positions. For those cases, bring 2 peptide templates, one for each position and combine them. 
- Change anchor prediction step based on the analysis done with the crude 3D rectangle boundaries
- Change anchor tolerance step to have options for Major_anchors/Secondary_anchors
- See to the random error of not writing files and not being able to find them
- Prepare the code for RMSD comparison with PANDORA paper.
- Testing/Testing/Testing...

### PTMs:

#### Major issues:
- Implement other PTMs:
	- Already done:
		- Phosphorylation
		- Citrullination
		- S-Nitrosylation
		- Methylation
		- Acetylation
		- Carbamylation
		- P-hydroxylation (non-biological 4S configuration not implemented)
		- C-Oxidation (not the R variation in CSX)
		- M-Oxidation (just the default variation)
	- To Do:
		- Nitration

#### Minor issues:
- Some PTMs are possible for any amino-acid in the N-termini. This has been observed, so future work can focus on adding this extra option. 
- Malondialdehyde adducts is particularly complex, so I'll leave it for now. 
- Deamindation is not needed for now, but since it has been observed in pMHCs a lot, maybe implement it by expanding `pytms.py`
- Some renaming on the extra hydrogens added by `pymol` will be necessary if we are to expand every PTM on openMM/GROMACS. 

### OpenMM:
- Update forcefield parameters for PTMs with the one released this year from Mackerel's group
	- https://www.charmm.org/charmm/resources/charmm-force-fields/
	- https://github.com/openmm/openmmforcefields/tree/master/charmm#Converting

This will probably happen with the new `openmmforcefields` release (when that happens). 

### Actual Experiments:

- RMSD comparison with PANDORA paper. Focus on peptides that are difficult to model in the database, as they don't have an obvious homologue. 
- Test RCD parameters/rigid vs. flexible/OpenMM vs. not OpenMM/#Ape-Gen cycles/any other experiment you can think of
- Performance on PTMed peptides

### Expansions:

- GUI
- Ensemble of receptor conformations in RCD step (maybe after that too?)
- Extend to Class-II using the same principle
- AlphaFold / RosettaTFold (wait for an API like thing) (could I also generate a bunch of HLAs offline and assess instead?)
- Think about alternative scoring function? )and test it's RMSD to vindardo/vina/CNN)
