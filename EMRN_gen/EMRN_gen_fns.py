import pandas as pd
import numpy as np
from tqdm import tqdm
import ast

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit import DataStructs
import re

# Function to canonicalize SMILES strings
def convert_to_canonical_smiles(smiles):
    if Chem.MolFromSmiles(smiles) is not None:
        cans = Chem.MolToSmiles(Chem.MolFromSmiles(smiles),True)
        return cans
    else:
        return smiles

def extract_uspto_reactants(id_reac):
    # id_reac['reactant list'] contains the reactants in the USPTO reactions
    reactant_strings = id_reac['reactant list']
    
    reactant_lists = reactant_strings.apply(ast.literal_eval)
    
    flattened_reactants = [item for sublist in reactant_lists for item in sublist]
    
    # List of the total indices: 0 - 1808937
    # id_reac['Unnamed: 0'] has theindices of the USPTO reactions as per the original USPTO database before filtering
    idx_list = list(id_reac['Unnamed: 0'])
    # Count of number of reactants per reaction in uspto db
    reac_count = list(reactant_lists.apply(lambda x: len(x)))
    
    # Based on reac_count, creating a list which has the index placed equal to the count for that index
    # Ex: idx_list = [0, 1, 2] and reac_count = [2, 3, 2]
    # then, full_idx_list = [0, 0, 1, 1, 1, 2, 2]
    full_idx_list = [elem for elem, count in zip(idx_list, reac_count) for _ in range(count)]
    
    data_dict = {}
    for elem, num in zip(flattened_reactants, full_idx_list):
        if elem not in data_dict:
            data_dict[elem] = [num]
        else:
            data_dict[elem].append(num)
    
    # Dataframe with all the unique reactants in id_reac with respective reaction indices they are present as reactants in
    unique_id_reac_df = pd.DataFrame(list(data_dict.items()), columns=['Reactant', 'Indices'])
    reac_series = unique_id_reac_df['Reactant']
    
    # Canonicalize the reactants in reac_series
    reac_series_can = pd.Series([convert_to_canonical_smiles(i) for i in tqdm(reac_series.values, desc='Processing')])

    return unique_id_reac_df, reac_series_can

# Function for USPTO reactant and metabolite-based precursor matching for a given round
def emrn_reactant_matching(met_smi, reac_series_can, unique_id_reac_df, id_reac):

    # Store the unique metabolite based reactants for a given round
    unique_round_reactants = pd.Series(met_smi.unique()).apply(convert_to_canonical_smiles)
    matching_smiles = reac_series_can[reac_series_can.isin(unique_round_reactants)]
    unq_smi_sub = unique_id_reac_df.loc[matching_smiles.index]

    id_list = []
    r_list = []
    p_list = []

    # Extracting the products for each of the matched reactants and storing all the info in lists
    for reac in range(len(unq_smi_sub)):
        x = id_reac.loc[id_reac['Unnamed: 0'].isin(unq_smi_sub['Indices'].iloc[reac])]
        prod_strings = x['product list']
        prod_lists = prod_strings.apply(ast.literal_eval)
        
        for i in range(len(prod_lists)):
            id_list.extend([x.iloc[i]['Unnamed: 0']]*len(prod_lists.values[i]))
            r_list.extend([unq_smi_sub.iloc[reac]['Reactant']]*len(prod_lists.values[i]))
            p_list.extend(prod_lists.values[i])

    # Dataframe to store reactant, product, reaction index
    # We write the Metabolite-derived reactant column as 'Met-reactant' in the rounds database
    r_net_df = pd.DataFrame({'Met_reactant': r_list, 'Product': p_list, 'Rxn_idx': id_list})
    # Removing reactions which have same molecule on reactant and product side
    r_net_df_filt = r_net_df.loc[~r_net_df['Rxn_idx'].isin(r_net_df.loc[r_net_df['Met_reactant']==r_net_df['Product']]['Rxn_idx'])]
    r_net_df_filt = r_net_df_filt.drop_duplicates().reset_index(drop=True)

    return r_net_df_filt

# Function for filtering by MW for 'sustainable reactions'
def mw_sus_filtering(r_net_df_filt, id_reac):

    r_smi = pd.Series(r_net_df_filt['Product'].unique())
    r = r_net_df_filt.loc[r_net_df_filt['Rxn_idx'].isin(id_reac['Unnamed: 0'])]
    r_source = r['Met_reactant']
    r_idx = r['Rxn_idx']
    
    r_source_can = [convert_to_canonical_smiles(i) for i in tqdm(r_source, desc='Processing')]

    # Precompute the molecular weights of the reactants
    source_mws = {source: Chem.Descriptors.MolWt(Chem.MolFromSmiles(source)) for source in set(r_source_can)}
    
    selected_rxn = []
    
    for idx, source in tqdm(zip(r_idx, r_source_can), total=len(r_idx), desc='Processing'):
        reac_list = ast.literal_eval(id_reac.loc[id_reac['Unnamed: 0'] == idx]['reactant list'].values[0])

        # Canonicalize reac_list
        reac_list_can = [convert_to_canonical_smiles(i) for i in reac_list]
        # Remove the source from the reactant list
        mod_reac_list = [i for i in reac_list_can if i != source]
    
        # Calculate molecular weights of modified reactant list in one line
        mw_reac_list = [Chem.Descriptors.MolWt(Chem.MolFromSmiles(i)) for i in mod_reac_list]
    
        # Use the precomputed molecular weight of the source
        source_mw = source_mws[source]
    
        # Check if all elements are smaller than the source molecular weight
        if all(mw < source_mw for mw in mw_reac_list):
            # Efficiently append the index
            matching_index = r.index[(r['Rxn_idx'] == idx) & (r['Met_reactant'] == source)]
            if not matching_index.empty:
                selected_rxn.append(matching_index[0])

    sustainable_rxns = r.loc[selected_rxn]

    return sustainable_rxns

def add_non_duplicates(combo, combined_df, new_df, round_tag):
    # Merge and identify duplicates
    if combo =='reaction':
        non_duplicates = new_df[~new_df[['Met_reactant', 'Product', 'Rxn_idx']].apply(tuple, 1).isin(
            combined_df[['Met_reactant', 'Product', 'Rxn_idx']].apply(tuple, 1)
        )]
    elif combo == 'nonreaction':
        non_duplicate = new_df[~new_df['smiles'].isin(combined_df['smiles'])]
    
    # Append only non-duplicate rows with the correct 'Round' tag
    combined_df = pd.concat([combined_df, non_duplicates], ignore_index=True)
    
    return combined_df
