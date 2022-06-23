echo Starting up...

cp -f /data/Docker_image/replacement_files/prepare_ligand4.py /conda/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py

grep -rl 'simtk.openmm' /conda/ | xargs sed -i 's/simtk.openmm/openmm/g'