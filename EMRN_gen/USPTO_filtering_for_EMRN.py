import pandas as pd
from tqdm import tqdm
import ast

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit import DataStructs
import re

def count_heavy_atoms(smiles):
    """Function to count heavy atoms (non-hydrogen) in a molecule given its SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    else:
        return 0  # If the molecule is invalid, return 0 heavy atoms

id_reac = pd.read_csv('identified_reactants.csv')
#uspto = pd.read_csv('USPTO_FULL.csv')

# Saving the indices of reactions with only one product
single_prod_rxn = [idx for idx in tqdm((pd.Series(range(len(id_reac)))), desc='Processing') if len(ast.literal_eval(id_reac.iloc[idx]['product list'])) == 1]
# Saving the indices of reactions with more than one product
nonsingle_prod_rxn = [idx for idx in tqdm((pd.Series(range(len(id_reac)))), desc='Processing') if len(ast.literal_eval(id_reac.iloc[idx]['product list'])) != 1]

id_reac_nonsingle_prod = id_reac.loc[nonsingle_prod_rxn]
id_reac_single_prod = id_reac.loc[single_prod_rxn]

# For the nonsingle products reactions, checking for reactions where only one products has > 4 heavy atoms
# making it a possible 'main reactant' of the reaction
reaction_indices = []
for idx in tqdm(id_reac['Unnamed: 0'].loc[nonsingle_prod_rxn], desc='Processing'):
    product_list = ast.literal_eval(id_reac.iloc[idx]['product list'])
    heavy_atom_counts = [count_heavy_atoms(prod) for prod in product_list]
    # Check if n-1 products have <= 3 heavy atoms, and one product has >= 4 heavy atoms
    if (len([count for count in heavy_atom_counts if count <= 3]) == len(product_list) - 1 and len([count for count in heavy_atom_counts if count >= 4]) == 1):
        reaction_indices.append(idx)

# Pick the reations from the nonsingle_prod_rxn list that contain the indices for 'main product' reactions
prod_filt = id_reac.loc[nonsingle_prod_rxn].loc[reaction_indices]
# Combine single product db and heavy-atom filtered db
id_reac_final = pd.concat([id_reac_single_prod, prod_filt]).sort_values(['Unnamed: 0']).reset_index(drop=True)

# Filtering out reactions which have 'unmodified products' - same molecules on reactant and product side
unmod_srs = id_reac_final['unmodified products list'].apply(ast.literal_eval)
idx = [index for index, i in enumerate(unmod_srs) if len(i) != 0]
id_reac_final_filtered = id_reac_final.drop(idx, inplace=True)

# Save final USPTO database to be used for EMRN generation
id_reac_final_filtered.to_csv('id_reac_final.csv')