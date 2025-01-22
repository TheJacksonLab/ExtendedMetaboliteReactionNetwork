import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import seaborn as sns
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import AllChem as Chem
import copy
from itertools import chain
import re
import ast
from rdkit.Chem import Draw
from rdkit.Chem import rdBase
from rdkit import RDConfig
import os
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

omg_mon = pd.read_csv('OMG_monomers.csv')
omg_pol = pd.read_csv('OMG_polymers.csv')

metr_sus_df = pd.read_csv('metr_sus_df.csv')
# Since we do this analysis for 2 reaction rounds
metr_sus_df2 = metr_sus_df.loc[metr_sus_df['tag'].isin(['Metabolite', 'Round1', 'Round2'])]

# Removing molecules that do not have import OMG reactant functional groups
# Store mol object in mol list and smiles in monomers_bag given that the mol objects are not None
# SMART annotation - dictionary with functional groups and their respective SMARTS

# List of OMG compatible functional groups
sub_structure_dict = {
            'acetylene': '[CX2]#[CX2]',
            'di_acid_chloride': '[CX3](=O)[Cl]',
            'conjugated_di_bromide': '[c;R][Br]',
            'cyclic_carbonate': '[OX1]=[CX3;R]([OX2;R][C;R])[OX2;R][C;R]',
            'cyclic_ether': '[C;R][O;R]([C;R])',  # '[OX2;R]([CX2;R][C;R])[CX2;R][C;R]'
            'cyclic_olefin': '[CH1;R][CH1;R]=[CH1;R][CH1;R]',
            'cyclic_sulfide': '[C;R][S;R]([C;R])',
            'di_amine': '[NX3H2;!$(NC=O)]',
            'di_carboxylic_acid': '[CX3](=O)[OX2H]',
            'di_isocyanate': '[NX2]=[CX2]=[OX1]',
            'di_ol': '[C,c;!$(C=O)][OX2H1]',
            'hydroxy_carboxylic_acid_OH': '[!$(C=O)][OX2H1]',
            'hydroxy_carboxylic_acid_COOH': '[CX3](=O)[OX2H]',
            'lactam': '[NH1;R][C;R](=O)',
            'lactone': '[O;R][C;R](=O)',
            'terminal_diene': '[CX3H2]=[CX3H1]',
            'vinyl': '[CX3;!R]=[CX3]'
        }


df2_list = list(metr_sus_df2['smiles'].unique())

# This code is taken from the Open Macromolecular Genome github and modified

mol = []
monomers_bag=[]
for smiles in df2_list:
    mol_obj = Chem.MolFromSmiles(smiles)
    if mol_obj is not None:
        mol.append(mol_obj)
        monomers_bag.append(smiles)
       
# find functional group components
df_sub_structure_unchanged = pd.DataFrame({'smiles': monomers_bag, 'mol': mol})
for key, value in sub_structure_dict.items():
    sub_structure_mol = Chem.MolFromSmarts(value)
    df_sub_structure_unchanged['%s' % key] = df_sub_structure_unchanged['mol'].apply(lambda x: len(x.GetSubstructMatches(sub_structure_mol)))
    
df_sub_structure = df_sub_structure_unchanged.copy()

# hydroxy_carboxylic_acid should have both OH and COOH
target = ['hydroxy_carboxylic_acid_OH', 'hydroxy_carboxylic_acid_COOH']
df_sub_structure['hydroxy_carboxylic_acid'] = df_sub_structure[target[0]] * df_sub_structure[target[1]]
df_sub_structure = df_sub_structure.drop(labels=[target[0], target[1]], axis=1)
# deal with exception
elim_idx = []
for col in df_sub_structure.columns:
    for row in df_sub_structure.index.values:
        # reduce di_* value to 0 if the value is 1. di_* value should be 2
        if 'di' in col and df_sub_structure.loc[row, col] == 1:
            df_sub_structure.loc[row, col] = 0
        # reduce di_* value to 0 if the value is larger than 2.
        if 'di' in col and df_sub_structure.loc[row, col] >= 3:
            df_sub_structure.loc[row, col] = 0
        # hydroxy_carboxylic_acid should have only one pair of OH and COOH
        if col == 'hydroxy_carboxylic_acid' and df_sub_structure.loc[row, col] >= 2:
            #print("[POLYMER] There is more than one hydroxy_carboxylic group in a monomer", flush=True)
            #return
            elim_idx.append(row)
        # if there are both vinyl group and cyclic olefin -> cyclic olefin has a priority (arbitrary)
        if col == 'vinyl' and df_sub_structure.loc[row, col] == 1 and \
                df_sub_structure.loc[row, 'cyclic_olefin'] == 1:
            df_sub_structure.loc[row, col] = 0
        # if there are both vinyl group and terminal diene -> terminal diene has a priority (arbitrary)
        if col == 'vinyl' and df_sub_structure.loc[row, col] == 2 and \
                df_sub_structure.loc[row, 'terminal_diene'] == 2:
            df_sub_structure.loc[row, col] = 0
        # number of vinyl functional group should be 1
        if col == 'vinyl' and df_sub_structure.loc[row, col] >= 2:
            #print("[POLYMER] There is more than one vinyl group in a monomer", flush=True)
            #return
            elim_idx.append(row)
        # number of acetylene functional group should be 1
        if col == 'acetylene' and df_sub_structure.loc[row, col] >= 2:
            #print("[POLYMER] There is more than one acetylene group in a monomer", flush=True)
            #return
            elim_idx.append(row)
        # if there are both cyclic ether and lactone -> lactone has a priority (arbitrary)
        if col == 'cyclic_ether' and df_sub_structure.loc[row, col] == 1 and \
            df_sub_structure.loc[row, 'lactone'] == 1:
            df_sub_structure.loc[row, col] = 0
        # if there are both cyclic ether and cyclic_carbonate -> cyclic carbonate has a priority (arbitrary)
        if col == 'cyclic_ether' and df_sub_structure.loc[row, col] == 2 and \
            df_sub_structure.loc[row, 'cyclic_carbonate'] == 1:
            df_sub_structure.loc[row, col] = 0
        # if there are both lactone and cyclic_carbonate -> cyclic carbonate has a priority (arbitrary)
        if col == 'lactone' and df_sub_structure.loc[row, col] == 2 and \
            df_sub_structure.loc[row, 'cyclic_carbonate'] == 1:
            df_sub_structure.loc[row, col] = 0
        # count the number of atoms in a ring of cyclic ether and cyclic sulfide
        # -> only 3, 4, and larger than 6 are allowed
        if (col == 'cyclic_ether' or col == 'cyclic_sulfide') and df_sub_structure.loc[row, col] >= 1:
            cyclic_smiles = df_sub_structure.loc[row, 'mol']
            ring_info = cyclic_smiles.GetRingInfo()
            for ring_cluster in ring_info.AtomRings():
                ring_idx = ring_cluster[0]
                min_num_atoms_in_ring = ring_info.MinAtomRingSize(ring_idx)
                if min_num_atoms_in_ring in (5, 6):
                    df_sub_structure.loc[row, col] = 0
        # number of functional group should be 1
        if 'di' not in col and col != 'smiles' and col != 'mol' and df_sub_structure.loc[row, col] >= 2:
            #print("[POLYMER] There is more than one cyclic functional group in a monomer", flush=True)
            #return
            elim_idx.append(row)

# decide if there are more than one functional groups in a molecule
reaction_sites = dict()
reaction_groups = dict()
reaction_monomers = dict()
for row_idx in range(df_sub_structure.shape[0]):
    count = 0
    for col in df_sub_structure.columns[2:]:  # exclude 'smiles' and 'mol' columns
        # not self-condensation of hydrocarboxylic acid
        if col != 'hydroxy_carboxylic_acid' and df_sub_structure.iloc[row_idx][col] != 0:
            count += 1
            # find react sites and reaction groups in SMILES molecules
            reaction_groups['monomer_%d' % (row_idx + 1)] = col
            reaction_monomers['monomer_%d' % (row_idx + 1)] = df_sub_structure.iloc[row_idx]['smiles']
            reaction_sites['monomer_%d' % (row_idx + 1)] = df_sub_structure.iloc[row_idx]['mol']. \
                GetSubstructMatches(Chem.MolFromSmarts(sub_structure_dict[col]))

        # self-condensation of hydroxycarboxylic acid
        elif col == 'hydroxy_carboxylic_acid' and df_sub_structure.iloc[row_idx][col] != 0:
            count += 1
            reaction_groups['monomer_%d' % (row_idx + 1)] = col
            reaction_monomers['monomer_%d' % (row_idx + 1)] = df_sub_structure.iloc[row_idx]['smiles']
            # add reaction sites of self-condensation of hyroxycarboxylic acid
            sites = ()
            for key, value in sub_structure_dict.items():
                if key.startswith(col):
                    sites += df_sub_structure.iloc[row_idx]['mol'].\
                        GetSubstructMatches(Chem.MolFromSmarts(value))
                    reaction_sites['monomer_%d' % (row_idx + 1)] = sites

    # check if there are more than two functional groups
    if count >= 2:
        #print("[POLYMER] There is more than one functional group in a monomer", flush=True)
        #return
        elim_idx.append(row)
    # check if there is no functional groups
    elif count == 0:
        #print("[POLYMER] There is a monomer that doesn't have a functional group", flush=True)
        #return
        elim_idx.append(row)

idx_to_eliminate = list(set(elim_idx))
# This df contains all the molecules in EMRN that have atleast one of the 17 OMG compatible functional groups
df_sub_structure_single = df_sub_structure.drop(idx_to_eliminate).reset_index(drop=True)

unq_pol_reac = pd.concat([omg_pol['reactant_1'], omg_pol['reactant_2']]).unique()
pol_overlap = metr_sus_df2.loc[metr_sus_df2['smiles'].isin(unq_pol_reac)]

# Finding exact matches in fg_df_single with the omg polymer reactants
omg_exact_match = df_sub_structure_single.loc[df_sub_structure_single['smiles'].isin(pol_overlap['smiles'])]
other_mols = df_sub_structure_single.loc[~df_sub_structure_single['smiles'].isin(pol_overlap['smiles'])]

# Filter the rows in metr_sus_df where 'smiles' is in df_sub_structure_single['smiles']
df_sub_structure_tag = metr_sus_df2.loc[metr_sus_df2['smiles'].isin(df_sub_structure_single['smiles']), 'tag'].tolist()
df_sub_structure_single['tag'] = df_sub_structure_tag