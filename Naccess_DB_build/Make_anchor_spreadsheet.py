import pandas as pd
from biopandas.pdb import PandasPdb
import numpy as np
import math

from tqdm import tqdm
import os
import glob

# z-dist of CA atoms to the beta sheet
anchor_list = []
for template_file in tqdm(os.listdir('./new_templates/pMHCI_AC/')):
    
    # Load template file
    template_file_name = './new_templates/pMHCI_AC/' + template_file
    template_peptide = PandasPdb()
    template_peptide.read_pdb(template_file_name)
    
    # Extract the CA atom locations of the peptide
    peptide = template_peptide.df['ATOM'][template_peptide.df['ATOM']['chain_id'] == 'C']
    peptide_length = max(peptide['residue_number'])
    peptide = peptide[['residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']]
    peptide_atom_list = list(peptide.itertuples(index=False, name=None))
    
    # Extract the CA atom locations of the beta sheet
    MHC = template_peptide.df['ATOM'][template_peptide.df['ATOM']['chain_id'] == 'A']
    beta_sheet = list(range(1, 45)) + list(range(95, 121))
    MHC = MHC[MHC['residue_number'].isin(beta_sheet)]
    MHC = MHC[MHC['atom_name'] == 'CA'][['residue_number', 'residue_name', 
                                         'x_coord', 'y_coord', 'z_coord']]
    MHC_atom_list = list(MHC.itertuples(index=False, name=None))
    
    # Initialization
    peptide_min_dist = {}
    for pos in range(1, peptide_length + 1):
        peptide_min_dist[pos] = float("inf")
    
    # Calculate minimum distance to the beta-sheet
    for peptide in peptide_atom_list:
        for MHC in MHC_atom_list:
            dist = math.sqrt((peptide[2] - MHC[2])**2 + (peptide[3] - MHC[3])**2 + (peptide[4] - MHC[4])**2)
            if dist < peptide_min_dist[peptide[0]]:
                peptide_min_dist[peptide[0]] = dist 
   
    # Major anchors:
    anchor_1 = np.argmin(list(peptide_min_dist.values())[:3]) + 1
    anchor_2 = np.argmin(list(peptide_min_dist.values())[-3:]) + peptide_length - 3 + 1
    major_anchors = str(anchor_1) + "," + str(anchor_2)
    anchor_list.append((template_file, anchor_1, anchor_2, peptide_length))

Beta_sheet_anchors_full = pd.DataFrame(anchor_list, columns=['pdb_code', 'beta_anchor_1_full', 'beta_anchor_2_full', 'peptide_length'])

# z-dist of all atoms (practically side_chains) to the beta sheet
anchor_list = []
for template_file in tqdm(os.listdir('./new_templates/pMHCI_AC/')):
    
    # Load template file
    template_file_name = './new_templates/pMHCI_AC/' + template_file
    template_peptide = PandasPdb()
    template_peptide.read_pdb(template_file_name)
    #peptide_sequence = template_peptide.amino3to1()
    #peptide_sequence = ''.join(peptide_sequence.loc[peptide_sequence['chain_id'] == 'C', 'residue_name'])
    
    # Extract the CA atom locations of the peptide
    peptide = template_peptide.df['ATOM'][template_peptide.df['ATOM']['chain_id'] == 'C']
    peptide_length = max(peptide['residue_number'])
    peptide = peptide[peptide['atom_name'] == 'CA'][['residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']]
    peptide_atom_list = list(peptide.itertuples(index=False, name=None))
    
    # Extract the CA atom locations of the beta sheet
    MHC = template_peptide.df['ATOM'][template_peptide.df['ATOM']['chain_id'] == 'A']
    beta_sheet = list(range(1, 45)) + list(range(95, 121))
    MHC = MHC[MHC['residue_number'].isin(beta_sheet)]
    MHC = MHC[MHC['atom_name'] == 'CA'][['residue_number', 'residue_name', 
                                         'x_coord', 'y_coord', 'z_coord']]
    MHC_atom_list = list(MHC.itertuples(index=False, name=None))
    
    # Initialization
    peptide_min_dist = {}
    for pos in range(1, peptide_length + 1):
        peptide_min_dist[pos] = float("inf")
    
    # Calculate minimum distance to the beta-sheet
    for peptide in peptide_atom_list:
        for MHC in MHC_atom_list:
            dist = math.sqrt((peptide[2] - MHC[2])**2 + (peptide[3] - MHC[3])**2 + (peptide[4] - MHC[4])**2)
            if dist < peptide_min_dist[peptide[0]]:
                peptide_min_dist[peptide[0]] = dist 
   
    # Major anchors:
    anchor_1 = np.argmin(list(peptide_min_dist.values())[:3]) + 1
    anchor_2 = np.argmin(list(peptide_min_dist.values())[-3:]) + peptide_length - 3 + 1
    major_anchors = str(anchor_1) + "," + str(anchor_2)
    anchor_list.append((template_file, anchor_1, anchor_2, peptide_length))

Beta_sheet_anchors_CA = pd.DataFrame(anchor_list, columns=['pdb_code', 'beta_anchor_1_CA', 'beta_anchor_2_CA', 'peptide_length'])

# If both z-dists indicate non-canonical anchor, it's all good, but if not, default to canonical
Beta_sheet_anchors = Beta_sheet_anchors_full.merge(Beta_sheet_anchors_CA, on=['pdb_code', 'peptide_length'])
Beta_sheet_anchors['beta_anchor_1'] = np.where(Beta_sheet_anchors['beta_anchor_1_full']==Beta_sheet_anchors['beta_anchor_1_CA'], Beta_sheet_anchors['beta_anchor_1_full'], 2)
Beta_sheet_anchors['beta_anchor_2'] = np.where(Beta_sheet_anchors['beta_anchor_2_full']==Beta_sheet_anchors['beta_anchor_2_CA'], Beta_sheet_anchors['beta_anchor_2_full'], Beta_sheet_anchors['peptide_length'])
Beta_sheet_anchors = Beta_sheet_anchors[['pdb_code', 'beta_anchor_1', 'beta_anchor_2', 'peptide_length']]

# RSA result
# An anchor is a major anchor if it has the lowest RSA
RSA = pd.read_csv("./RSA_csvs/RSA.csv", header=None)
RSA.columns = ['pdb_code', 'AA', 'AA_index', 'SASA_All', 'RSA_all', 'SASA_sidechain', 'RSA_sidechain']
RSA_anchor_1 = RSA.loc[RSA.groupby("pdb_code").head(3).groupby(['pdb_code'])["RSA_sidechain"].idxmin()][['pdb_code', 'AA_index']]
RSA_anchor_1.columns = ['pdb_code', 'RSA_anchor_1']
RSA_anchor_2 = RSA.loc[RSA.groupby("pdb_code").tail(3).groupby(['pdb_code'])["RSA_sidechain"].idxmin()][['pdb_code', 'AA_index']]
RSA_anchor_2.columns = ['pdb_code', 'RSA_anchor_2']
RSA_major_anchors = RSA_anchor_1.merge(RSA_anchor_2, on=['pdb_code'])

# Make a consensus of these approaches
Consensus_anchors = Beta_sheet_anchors.merge(RSA_major_anchors, on=['pdb_code'])
Consensus_anchors['consensus_anchor_1'] = np.where(Consensus_anchors['beta_anchor_1']==Consensus_anchors['RSA_anchor_1'], Consensus_anchors['beta_anchor_1'], '2')
Consensus_anchors['consensus_anchor_2'] = np.where(Consensus_anchors['beta_anchor_2']==Consensus_anchors['RSA_anchor_2'], Consensus_anchors['beta_anchor_2'], Consensus_anchors['peptide_length'])
Consensus_anchors['Major_anchors'] = Consensus_anchors['consensus_anchor_1'].astype(str)+','+Consensus_anchors['consensus_anchor_2'].astype(str)
Consensus_anchors = Consensus_anchors[['pdb_code', 'Major_anchors']]

# Calculate also Secondary anchors given thresholds in literature:
RSA = RSA[RSA['RSA_sidechain'] <= 20]
RSA = RSA[((RSA['AA'].isin(['ALA', 'GLY', 'SER'])) & (RSA['RSA_all'] <= 30)) | (~RSA['AA'].isin(['ALA', 'GLY', 'SER']))]
RSA = RSA.rename(columns={"AA_index": "Secondary_anchors"})
RSA['Secondary_anchors'] = RSA['Secondary_anchors'].astype(str)
RSA_Secondary_anchors = RSA.groupby(['pdb_code'])['Secondary_anchors'].apply(','.join)
Consensus_anchors = Consensus_anchors.merge(RSA_Secondary_anchors, on=['pdb_code'])

# Unify Major and Secondary anchors if some Major anchor got away
Major_anchors = Consensus_anchors['Major_anchors'].tolist()
Major_anchors = [set(x.split(",")) for x in Major_anchors]
Secondary_anchors = Consensus_anchors['Secondary_anchors'].tolist()
Secondary_anchors = [set(x.split(",")) for x in Secondary_anchors]
for i, x in enumerate(Secondary_anchors):
    Secondary_anchors[i] = Secondary_anchors[i].union(Major_anchors[i])
Secondary_anchors = [",".join(x) for x in Secondary_anchors]
Consensus_anchors['Secondary_anchors'] = Secondary_anchors

# Add pdb_code and allele allotype
pdb_code_list = []
hla_template_list = []
with open("./template_sequences/template_names.log", "r") as f:
    while True:
        line1 = f.readline().strip("\n")
        line2 = f.readline().strip("\n")
        if not line2: 
            break  # EOF
        else:
            pdb_code_list.append(line1)
            hla_template_list.append(line2)
d = {'pdb_code': pdb_code_list, 'MHC': hla_template_list}
df = pd.DataFrame(data=d)

# Add peptide sequence
pdb_code_list = []
peptide_sequence_list = []
for filename in glob.glob('./template_sequences/peptide/*.pdb'):
    with open(filename, "r") as f:
        while True:
            line1 = f.readline().strip("\n").replace(">", "")
            line2 = f.readline().strip("\n")
            if not line2: 
                break  # EOF
            else:
                pdb_code_list.append(line1)
                peptide_sequence_list.append(line2)
d = {'pdb_code': pdb_code_list, 'peptide': peptide_sequence_list}
df_peptide = pd.DataFrame(data=d)
df_peptide['peptide_length'] = df_peptide['peptide'].apply(len)
df = df.merge(df_peptide, on = ["pdb_code"])


Template_information = df.merge(Consensus_anchors, on=['pdb_code'])
Template_information.to_csv("./Template_DB_information.csv")