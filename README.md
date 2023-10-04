# New Ape-Gen version:

## Prelims

To set up APE-Gen2.0 locally, download the repo using `git clone`.

Before you move on to anything else, am important step must be performed: Download the frequency matrices generated using MHCFlurry2.0 through the link below:

https://rice.box.com/s/duxshqxtkykg7u9y3b0tl79j8byc2dp2

These will be crucial for the softawre to run. After these are downloaded, go ahead and put those on the `helper_files` folder. Then you're all set for the installation!

## Installation

We recommend that APE-Gen2.0 is set up with a docker configuration, so the latest version of `docker` should be installed on your system. Setting APE-Gen2.0 through a python virtual env is definitely plausible (check the `environment.yml`and the provided `Dockerfile` for all needed packages), but we'll leave that to you if you want to set it up this way. 

Just running:
```
bash build_image.sh
```

should build the docker image After this, just running:
```
bash start_image.sh
```

gets you into the container, APE-Gen2.0 can be run from there. 

## Usage

To run APE-Gen2.0 for existing MHC templates in the database, you should be good to go (Ignore the warning message about the MODELLER invalid key). However, for a new MHC not found in the database, you will need to provide a valid MODELLER key. To do that, go to:

```
cd /conda/envs/apegen/lib/modeller-10.4/modlib/modeller/
```

and modify the lines of `config.py` with a valid MODELLER key.

For now, peptide sequence (along with positions with PTM) + an MHC allotype should work! (more inputs to follow...)

```
python New_APE-Gen.py LLGIGSLTV HLA-A*02:01 --verbose
python New_APE-Gen.py LLGIGSLTV HLA-A*02:01 --verbose --score_with_openmm
python New_APE-Gen.py LLGIGpSLTV HLA-A*02:01 --verbose --score_with_openmm
python New_APE-Gen.py LLGpSGpSLTV HLA-A*02:01 --verbose --score_with_openmm
```

Please contact the team should you have any issues with running APE-Gen2.0!