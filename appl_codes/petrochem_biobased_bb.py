import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit import DataStructs

def convert_to_canonical_smiles(smiles):
    if Chem.MolFromSmiles(smiles) is not None:
        cans = Chem.MolToSmiles(Chem.MolFromSmiles(smiles),True)
        return cans
    else:
        return smiles

# Function to get parts of the EMRN dataframe that contain a certain list subset of building block chemicals
def get_bb_df(bb_list, met_df):
    bb_list_can = [convert_to_canonical_smiles(i) for i in bb_list]
    final_df = met_df.loc[met_df['smiles'].isin(bb_list_can)]
    return final_df

# Different building block subsets
# Petrochemical basic building blocks
pet_bb = ['C=C', 'CC=C', 'C=CCC', 'C=CC=C', 'CC(=C)C', 'C', 'C1=CC=CC=C1',
         'CC1=CC=CC=C1', 'CC1=CC=CC=C1C', 'CC1=CC(=CC=C1)C', 'CC1=CC=C(C=C1)C']

# Derivatives of basic petrochemical building blocks
pet_der = ['CCO', 'C1CO1', 'CC(=O)OC=C', 'CC=O', 'C(CCl)Cl', 'C=CC(=O)O', 'C=CC#N', 'CC(C)O', 'CC1CO1', 'CCCC=O', 'C=CCCl', 'CO', 'CCC1=CC=CC=C1',
          'CC(C)C1=CC=CC=C1', 'C1CCCCC1', 'C1=CC=C(C=C1)N+[O-]', 'CC1=C(C=C(C=C1)N=C=O)N=C=O', 'C1=CC=C(C=C1)C(=O)O']

# Bio-based platform chemicals
doe_pc = ['C(=CC(=O)O)C(=O)O','C(C(C(=O)O)O)C(=O)O' ,'C1=C(OC(=C1)C(=O)O)C(=O)O' ,'C(CO)C(=O)O','C(C(C(=O)O)N)C(=O)O','C(C(C(C(=O)O)O)O)(C(C(=O)O)O)O','C(CC(=O)O)C(C(=O)O)N' ,'C=C(CC(=O)O)C(=O)O','CC(=O)CCC(=O)O','C1C(COC1=O)O','C(C(CO)O)O','C(C(C(C(C(CO)O)O)O)O)O','C(C(C(C(CO)O)O)O)O','C(CC(=O)O)C(=O)O', 'OCC', 'C1=COC(=C1)C=O', 'c1cc(oc1CO)C=O', 'CC(=C)C=C', 'CC(C(=O)O)O', 'C(CO)C=O']

# Popular commodity monomers
com_mon = ['C=C', 'CC=C', 'C=CC1=CC=CC=C1', 'C=CCl', 'C1=CC(=CC=C1C(=O)O)C(=O)O', 'C(CO)O', 'C=CC#N', 'C=CC=C', 'C1CCC(=O)NCC1']

# met_df = pd.read_csv('metr_sus_df.csv')

