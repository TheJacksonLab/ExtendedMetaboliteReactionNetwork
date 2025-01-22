import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit import DataStructs
import scscore as scs

# Example of calculating SCScore and molecular weight for a set of SMILES
met = pd.read_csv('met_all_noMW_common.csv')
# SCScore calculation:

# load SC score class
model = scs.SCScorer()
model.restore()

# calculate SC Score
met['SC_score'] = list(met['canonical smiles'].apply(lambda x: model.get_score_from_smi(x)[1]))

# Molecular weight calculation:
mw_reac_list = [Chem.Descriptors.MolWt(Chem.MolFromSmiles(i)) for i in met['canonical smiles']]
met['Mol_Wt'] = mw_reac_list