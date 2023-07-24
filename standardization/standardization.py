#install rdkit
pip install rdkit


#importing the cheminformatics module needed
import molvs
import pandas as pd
import numpy as np
from molvs import Standardizer, normalize
from rdkit import Chem
from tqdm import tqdm
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.drawOptions.comicMode=True
import rdkit
print(rdkit.__version__)

#load the datasets saved using dataframe.to_csv('Replicase_bioactivity_data_preprocessed.csv', index=False)
df =pd.read_csv("C:/Users/ogbue/Replicase_bioactivity_data_preprocessed.csv")

#Define a canonical function to canonicalize the smiles
def Canonical(smiles):
  canon_smi = [Chem.CanonSmiles(smi) for smi in smiles]
  return canon_smi
# Generate the canonical smiles
canon_smiles = Canonical(df.canonical_smiles)
#replace the smiles column with canon smiles
df['canonical_smiles'] = canon_smiles
df
df.shape


# Create a Standardizer object
standardizer = Standardizer()

# Modify the normalization rules
norms = list(normalize.NORMALIZATIONS)
for i in range(len(norms) - 1, 0, -1):
    if norms[i].name == "Sulfoxide to -S+(O-)":
        del norms[i]
norms.append(normalize.Normalization("[S+]-[O-] to S=O", "[S+:1]([O-:2])>>[S+0:1](=[O-0:2])"))

# Set the modified normalization rules
standardizer.normalizations = norms

# Set the prefer_organic option to True
standardizer.prefer_organic = True
#standardize the canonical smiles and the rdkit mol object
stand_mol = []
stand_smi = []
for smi in df['canonical_smiles'].tolist():
    mol = Chem.MolFromSmiles(smi)  # Convert SMILES to RDKit Mol object
    try:
        standardized_mol = standardizer.standardize(mol)
        stand_mol.append(standardized_mol)
        standardized_smi = Chem.MolToSmiles(standardized_mol)
        stand_smi.append(standardized_smi)
        print(stand_smi)
    except Exception as e:
        print("Validation error: ", e)
            
#put the standardized smiles in the original dataframe
df['standardized_smiles'] = pd.DataFrame(stand_smi)
df.shape
#reset the index
df.reset_index(drop=True, inplace=True)
df.to_csv("Replicase_stand_smi_data.csv", index=False)
