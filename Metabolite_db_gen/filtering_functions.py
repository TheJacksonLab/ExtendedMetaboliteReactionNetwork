import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Function to canonicalize SMILES strings
def convert_to_canonical_smiles(smiles):
    if Chem.MolFromSmiles(smiles) is not None:
        cans = Chem.MolToSmiles(Chem.MolFromSmiles(smiles),True)
        return cans
    else:
        return smiles

# FILTERING
# Functions for data filtering for metabolites

# (i) only contain atoms frequently appearing in organic compounds (e.g. H, B, C, N, O, F, Al, Si, P, S, Cl, and Br), <br>
# (ii) not contain isotopes (e.g. 13C), <br>
# (iii) not contain explicit hydrogens (e.g. [H]OCCO),<br>
# (iv) not contain non-zero formal charges (e.g. COC(=O)c1ccccc1[Br+]c1ccccc1), <br>
# (v) not contain a period punctuation mark to include more than one molecule (e.g. Br.Oc1cc(on1)C1CCNCC1), and <br>
# (vi) not contain atom-map indices (e.g. [H:4][n:3]1[cH:2][cH:1][c:6]([CH2:7][Cl:8])[n:5]1).

# Limiting atomic number
def limit_atomic_num(mol):

    keep_list = [1, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 35]

    num_atm = mol.GetNumAtoms()
    flag = 0
    for i in range(num_atm):
        atomic_num = mol.GetAtomWithIdx(i).GetAtomicNum()
        if atomic_num not in keep_list:
            flag = 1
            break
    return flag

# Excluding isotopes
def limiting_isotope(mol):
    cnt = 0
    atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]
    for atom, isotope in atom_data:
        if isotope:
            cnt += 1
            break
    return cnt

# Metabolite database contains a column containing SMILES strings of the metabolites, with metabolite tag also.

def filtering_database(df, smiles_column):
    
    df['mol'] = df[smiles_column].apply(lambda x: Chem.MolFromSmiles(x))
    df = df.dropna().reset_index(drop=True)
    
    # Limiting atomic numbers
    df['flag'] = df['mol'].apply(lambda x: limit_atomic_num(x))
    df = df[df['flag'] == 0]
    
    # Checking for periods in smiles
    df['point'] = df[smiles_column].apply(lambda x: '.' in x)
    df = df[~df['point']]
    
    # Excluding [H] and atom mappings
    df['flag'] = df[smiles_column].apply(lambda x: ('[H]' in x) | (':' in x))
    df = df[df['flag'] == 0]
    
    # Exclude isotopes
    df['flag'] = df['mol'].apply(lambda x: limiting_isotope(x))
    df = df[df['flag'] == 0]
    
    # Excluding formal charges and radicals
    # Remove non-zero Formal charge
    df['flag'] = df['mol'].apply(lambda x: Chem.GetFormalCharge(x))
    df = df[df['flag'] == 0]
    
    # Removing '+' ions
    df['flag'] = df[smiles_column].apply(lambda x: '+' in x)
    df = df[df['flag'] == 0]
    
    # Removing '-' ions
    df['flag'] = df[smiles_column].apply(lambda x: '-' in x)
    df = df[df['flag'] == 0]
    
    # Removing radicals
    df['flag'] = df['mol'].apply(lambda x: Descriptors.NumRadicalElectrons(x))
    df = df[df['flag'] == 0]

    return df

