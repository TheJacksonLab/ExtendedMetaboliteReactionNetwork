import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import filtering_functions as ff

# Loading the E. coli, S. cerevisiae and P. aeruginosa databases
ecoli = pd.read_json('ecmdb.json')
pmdb = pd.read_excel('PaMet.xlsx')
yeast = pd.read_csv('yeast_detected_and_quantified.csv')

# Add a 'tag' column and combine the SMILES and tag columns from each database
ecoli_sm = ecoli.assign(tag='ecoli')[["moldb_smiles", "tag"]]
yeast_sm = yeast.assign(tag='yeast')[["SMILES", "tag"]]
pmdb_sm = pmdb.assign(tag='pmdb')[["SMILES", "tag"]]

# Rename columns for consistency and combine all datasets
met = (
    pd.concat([
        ecoli_sm.rename(columns={"moldb_smiles": "SMILES"}), 
        yeast_sm, 
        pmdb_sm
    ], axis=0, ignore_index=True)
    .drop_duplicates()
    .dropna()
    .reset_index(drop=True)
)

# Convert SMILES to canonical SMILES and create a new DataFrame
met_df = met.copy()
met_df["canonical smiles"] = met_df["SMILES"].apply(ff.convert_to_canonical_smiles)

# Retain only the required columns
met_df = met_df[["canonical smiles", "tag"]]

# Assign 'common' tag to rows with duplicated canonical smiles - these repesent the metabolites common between 2 or more of the organisms
met_df.loc[met_df["canonical smiles"].duplicated(keep=False), "tag"] = "common"

met_df = met_df.drop_duplicates().reset_index(drop=True)

# Carry out filtering of the metabolite database based on the following criteria:
# (i) only contain atoms frequently appearing in organic compounds (e.g. H, B, C, N, O, F, Al, Si, P, S, Cl, and Br), <br>
# (ii) not contain isotopes (e.g. 13C), <br>
# (iii) not contain explicit hydrogens (e.g. [H]OCCO),<br>
# (iv) not contain non-zero formal charges (e.g. COC(=O)c1ccccc1[Br+]c1ccccc1), <br>
# (v) not contain a period punctuation mark to include more than one molecule (e.g. Br.Oc1cc(on1)C1CCNCC1), and <br>
# (vi) not contain atom-map indices (e.g. [H:4][n:3]1[cH:2][cH:1][c:6]([CH2:7][Cl:8])[n:5]1).

met_df = ff.filtering_database(met_df, 'canonical smiles')

# Save as .csv
#met_df.to_csv("met_all_noMW_common.csv")