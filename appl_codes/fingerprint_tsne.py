import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import tqdm

import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, GraphDescriptors, Descriptors
from rdkit.Chem import Draw

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import openTSNE

# Calculating Morgan fingerprints
def morgan_fp(smi_list):
    smi_mols = [Chem.MolFromSmiles(x) for x in smi_list]
    
    radius = 3
    nBits = 1024 
    info = {}
    
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits, bitInfo=info) for mol in smi_mols]
    fps_np_array = []
    print('There are {} fingerprints'.format(len(info)))
    
    for fp_object in fps:
        fp_vect = np.zeros((1,), dtype=float)
        DataStructs.ConvertToNumpyArray(fp_object, fp_vect)
        fps_np_array.append(fp_vect)
        
    return fps_np_array

# Calculating tSNE coordinates
def tsne_calc(fps_array):
 
    pca_50 = PCA(n_components=50)
    crds_50 = pca_50.fit_transform(fps_array)
    
    %time crds_embedded = openTSNE.TSNE(n_components=2, perplexity=70, verbose=True, n_iter=3000).fit(crds_50)
    crds_tsne_df = pd.DataFrame(crds_embedded,columns=["X","Y"])

    return crds_tsne_df