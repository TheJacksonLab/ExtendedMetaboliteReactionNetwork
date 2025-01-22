import pandas as pd
import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, GraphDescriptors, Descriptors
from rdkit.Chem import Draw

# SMART annotation - dictionary with functional groups and their respective SMARTS - some popular functional groups in medicinal chem

sub_structure_dict = {
            'Alkene': '[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]',
            'Alkyne': '[CX2]#[CX2]',
            'Aldehyde': '[CX3H1](=O)[#6]',
            'Alkyl halide': '[#6;!R][F,Cl,Br,I]',
            'Aryl halide': '[#6;R][F,Cl,Br,I]',
            'Acyl halide': '[#6](=[O])[F,Cl,Br,I]',
            'Anhydride': '[CX3](=[OX1])[OX2][CX3](=[OX1])',
            'Amines': '[NX3+0,NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]',
            'Amide-C': '[NH2][CX3](=[OX1])[#6]',
            'Amide-N': '[NH2][CX3](=[OX1])[#7]',
            'Nitrile': '[NX1]#[CX2]',
            'Ketone': '[CX3;$(C([#6])(=[O])[#6])] (=[O;!$([O][O])])',
            'Ether': '[OD2]([#6;!R;!$(C=O)!$(CC)])[#6;!R;!$(C=O)]',
            'Alcohol': '[CX4][OH]',  # '[OX2;R]([CX2;R][C;R])[CX2;R][C;R]' # [OX2H][CX4;!$(C([OX2H])[O,S,#7,#15])]
            'Ester': '[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]',
            'Carboxylic acid': 'O=C[OH]',
            'Phosphate': '[PX4](=[OX1])([OX2H,OX1-,OX2R])([OX2H,OX1-,OX2R])([OX2H,OX1-,OX2R])'     
        }

# Checking for presence of above functional groups in a list of SMILES strings
def fg_mols(smi_list):
    # Store mol object in mol list and smiles in monomers_bag given that the mol objects are not None
    mol = []
    monomers_bag=[]
    for smiles in smi_list:
        mol_obj = Chem.MolFromSmiles(smiles)
        if mol_obj is not None:
            mol.append(mol_obj)
            monomers_bag.append(smiles)
            
    # find functional group components
    df_sub_structure = pd.DataFrame({'smiles': monomers_bag, 'mol': mol})
    for key, value in sub_structure_dict.items():
        sub_structure_mol = Chem.MolFromSmarts(value)
        df_sub_structure['%s' % key] = df_sub_structure['mol'].apply(lambda x: len(x.GetSubstructMatches(sub_structure_mol)))
    
    df_sub_structure['tag'] = met['tag']
    # Total count for each functional group
    #fg = list(df_sub_structure.iloc[:,2:-1].sum(axis=0))

    # Total count of number of molecules containing respective functional groups
    fg = list((df_sub_structure.iloc[:,2:-1] > 0).sum())
    
    # List of all functional groups checked for
    fg_cols = list(df_sub_structure.iloc[:,2:-1].columns)

    return fg, fg_cols