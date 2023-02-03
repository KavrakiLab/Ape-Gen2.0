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
cd /conda/envs/apegen/lib/modeller-10.4/modlib/modeller/
```

and modify the lines of `config.py` with a valid MODELLER key. Contact the lab to get a valid key.

Finally, download the files from this link:

```
https://rice.box.com/s/q4bh18eqvmvpzulqyfoszebo9jfkivc3
```

and put them in the helper files folder, as these are mandatory for the anchor prediction step!


For now, peptide sequence (along with phosphorylated positions) + HLA allotype works, but native + other types of inputs will follow:

```
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --verbose
python New_APE-Gen.py ARSEDEVILS HLA-A*11:01 --verbose --score_with_openmm
python New_APE-Gen.py ARpSEpTEVIpYS HLA-A*11:01 --verbose --score_with_openmm
```

## TO-DOs:

### Main Workflow:

#### Minor issues:
- Motifs from openvax that are identical in regards to alleles do not exist as is in the file. I need to match e.g. `A*02:17` to `A*02:01` correctly.  
- See if there is a cleaner way to define the anchors for cases not found in PMBEC.
- Add MHCFlurry motifs prediction in the workflow (you can take those from the notebooks created)
- Add pdb2pqr for protonation.
- Taxonomy on what is a macro and what is a class method must be done at some point
- Send an e-mail to Chacon lab for the RCD file that does not want to work
- Wrap the anchor extraction/other things that can be wrapped in the Peptide/pMHC class from the initial part up to the RCD (but also in general) and put them in the macro file, they look ugly.
- Put a check so that if a peptide template with an MHC that is the same as the input allotype is found, don't bother fetching smth different, and just assign the same template as the receptor (that would require flipping the template selection order though). 
- `_A`s and `_B`s need to go from the templates, I should be keeping only the As (rm files from history)
- I need to clean the dataframe fields. There are many that are basically useless, or I compute many while executing the algorithm, and pre-computation of those would be much better. 
- Convert macros to the way they previously were, or make them class methods if they concern exclusively a class. 
- `successful_conformations_statistics.csv` correction to `successfull`
- Eventually prune image space by removing unwanted packages from `environment.yml`
- Add all the folders with the executables to `$PATH`, so that we do not have absolute paths to Autodocktools etc. (if possible)
- Maybe have an argument for `autobox` size for SMINA scoring (larger values will result in greater time but much more accurate docking)
- Possibly insert GNINA in workflow somehow? (and test it's RMSD to vinardo/vina).
- Add `num_of_rounds` in homology modelling as an argument. 
- Make function in the places where there is repeated code
- Add `number_of_tries` parameter to OpenMM
- CSP routine that builds the correct atom names in the flexible file could probably be more optimized (CSP could potentially solve the whole thing instead of per residue?)
- The `receptor_template_file` variable in Section 1b is questionable, investigate more
- Consider peptide methods which receptor is an input to be transferred to the `pMHC_class`, and a pMHC object defined
- REVISE `simtk.openmm` import warning coming from PDBFixer, wait for a new version or modify the .xml: https://github.com/openmm/pdbfixer/issues/233
- Function documentation needs to be done thouroughly at some point
- Thorough input checking (example is peptide sequence in HLA peptide fetching must be an amino acid sequence)

#### Major issues:
- Testing/Testing/Testing...

### PTMs:

#### Major issues:
- Implement other PTMs (DONE!):
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
		- Nitration

#### Minor issues:
- Some PTMs are possible for any amino-acid in the N-termini. This has been observed, so future work can focus on adding this extra option. 
- Malondialdehyde adducts is particularly complex, but maybe could be added in the future. 
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
- Extend to Class-II using the same principle
- Think about alternative scoring function? (and test it's RMSD to vindardo/vina/CNN)
