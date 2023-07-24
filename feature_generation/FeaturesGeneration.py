# install rdkit
pip install rdkit
#install data progress bar library
pip install tqdm
#importing the cheminformatics module needed
import molvs
import pandas as pd
import numpy as np
from rdkit import Chem
import tempfile
import os
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem import MolStandardize
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from tqdm import tqdm
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.drawOptions.comicMode=True
import rdkit
print(rdkit.__version__)

#load the datasets containing the replicase standardized smiles
df = pd.read_csv('C:/Users/ogbue/Replicase_stand_smi_data.csv')
df.head()
#generate the standardized molecule from the standardized smiles 
stand_mol = []
stand_smi = df["standardized_smiles"].tolist()
for i in stand_smi:
    standardized_mol = Chem.MolFromSmiles(i)
    stand_mol.append(standardized_mol)
df['standardized_molecule'] = pd.DataFrame(stand_mol)
#generate the maccs keys chemical fingerprints for substructure search
def generate_MACCSfpts(data):
    maccs_fpts = []
    for mol in tqdm(data):
        mkeyfpts = MACCSkeys.GenMACCSKeys(mol)
        maccs_fpts.append(mkeyfpts)
    return np.array(maccs_fpts)
maccs_fpts = generate_MACCSfpts(df['standardized_molecule'])
maccskeys_fingerprints = pd.DataFrame(maccs_fpts, columns=['Col_A_{}'.format(i + 1)
                                  for i in range(maccs_fpts.shape[1])])
maccskeys_fingerprints

#path based fingerprints - rdkfingerprints, Daylight-like
def generate_RDKfpts(data):
    RDK_fpts = []
    for mol in tqdm(data):
        rdkfpts = AllChem.RDKFingerprint(mol, maxPath=5, fpSize=2048, nBitsPerHash=2 )
        RDK_fpts.append(rdkfpts)
    return np.array(RDK_fpts)
RDK_fpts = generate_RDKfpts(df['standardized_molecule'])
#put it in dataframe
RDK_fingerprints = pd.DataFrame(RDK_fpts, columns=['Col_B_{}'.format(i + 1)
                                  for i in range(RDK_fpts.shape[1])])
##toplogical fingerprints - atom pair
def generate_APfpts(data):
    AP_fpts = []
    for mol in tqdm(data):
        apfpts = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=2048)
        AP_fpts.append(apfpts)
    return np.array(AP_fpts)

AP_fpts = generate_APfpts(df['standardized_molecule'])
#put it in dataframe
AP_fingerprints = pd.DataFrame(AP_fpts, columns=['Col_C_{}'.format(i + 1)
                                  for i in range(AP_fpts.shape[1])])
AP_fingerprints.head()
AP_fingerprints.shape
##toplogical fingerprints - topological torsion
def generate_TTfpts(data):
    TT_fpts = []
    for mol in tqdm(data):
        ttfpts = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=2048)
        TT_fpts.append(ttfpts)
    return np.array(TT_fpts)

TT_fpts = generate_TTfpts(df['standardized_molecule'])
#put it in dataframe
TT_fingerprints = pd.DataFrame(TT_fpts, columns=['Col_D_{}'.format(i + 1)
                                  for i in range(TT_fpts.shape[1])])
TT_fingerprints.head()
TT_fingerprints.shape
##extended connectivity finger prints -  FINGERPRINTS
def generate_Morganfpts(data, radius):
    MORGAN_fpts = []
    for mol in tqdm(data):
        morganfpts = AllChem.GetMorganFingerprintAsBitVect(mol,radius, nBits=2048)
        MORGAN_fpts.append(morganfpts)
    return np.array(MORGAN_fpts)

ECFP4= generate_Morganfpts(df['standardized_molecule'], 2)
#put it in dataframe
ECFP4_fingerprints = pd.DataFrame(ECFP4, columns=['Col_E_{}'.format(i + 1)
                                  for i in range(ECFP4.shape[1])])
ECFP4_fingerprints.head()
ECFP4_fingerprints.shape
##extended connectivity finger prints - FEATURE CONNECTIVITY FINGERPRINTS
def generate_Featurefpts(data, radius):
    FEATURE_fpts = []
    for mol in tqdm(data):
        featurefpts = AllChem.GetMorganFingerprintAsBitVect(mol,radius, useFeatures=True, nBits=2048)
        FEATURE_fpts.append(featurefpts)
    return np.array(FEATURE_fpts)

FCFP4= generate_Morganfpts(df['standardized_molecule'], 2)
#put it in dataframe
FCFP4_fingerprints = pd.DataFrame(FCFP4, columns=['Col_F_{}'.format(i + 1)
                                  for i in range(FCFP4.shape[1])])
FCFP4_fingerprints.head()
FCFP4_fingerprints.shape

#calculate the RDKit descriptors
def RDKit_descriptors(smiles):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_names = calc.GetDescriptorNames()
    mol_descriptors = []
    for mol in mols:
        mol = Chem.AddHs(mol)
        descriptors = calc.CalcDescriptors(mol)
        mol_descriptors.append(descriptors)
    return mol_descriptors, desc_names
mol_descriptors, desc_names = RDKit_descriptors(stand_smi)
df_with_200descriptors = pd.DataFrame(mol_descriptors, columns=desc_names)
df_with_200descriptors 
#generate the toplological pharmacophore atom triplets fingerprints
#first convert the standardized smiles to an sdf file and generate file path
def toSDF(smiles_list):
    temp_dir = tempfile.mkdtemp()
    w = Chem.SDWriter(os.path.join(temp_dir, "temp.sdf"))

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            AllChem.Compute2DCoords(mol)
            w.write(mol)

    w.close()

    sdf_path = os.path.join(temp_dir, "temp.sdf")
    return sdf_path
#Generate the SDF file using the function toSDF
sdf_file_path = toSDF(stand_smi)
print("Generated SDF file path:", sdf_file_path)
#now copy the sdf_file path
#Download and install strawberry perl from https://strawberryperl.com/ to use perl script
#run the following command on the terminal to generate a csv file polyproteinTPATF.csv containing the tpatf.
# perl TopologicalPharmacophoreAtomTripletsFingerprints.pl --AtomTripletsSetSizeToUse FixedSize -v ValuesString -r polyproteinTPATF -o sdf_file_path

#open the csv file and put the fingerprints into array
file_path = #"file_path..to../polyproteinTPATF.csv"

with open(file_path, 'r') as f:
    all_features = []
    for line in f.readlines():
        if "Cmpd" in line:
            line = line.split(';')[5].replace('"', '')
            features = [int(i) for i in line.split(" ")]
            all_features.append(features)

features_array = np.array(features)
all_features_array = np.array(all_features)
all_features_array.shape

all_features_array

TPAT_fingerprints = pd.DataFrame(all_features_array, columns=['Col_G_{}'.format(i + 1)
                                  for i in range(all_features_array.shape[1])])
TPAT_fingerprints.head()
#TPAT_fingerprints.shape

# Concatenate the data frames column-wise
result = pd.concat([df, maccskeys_fingerprints,RDK_fingerprints, AP_fingerprints, TT_fingerprints, ECFP4_fingerprints, TPAT_fingerprints, df_with_208descriptors], axis=1)
result.shape

# Step 1: Check for NaN or blank spaces in all columns
rows_with_nan_or_blank = result.isnull().any(axis=1)

# Step 2: Drop the rows with NaN or blank spaces from the DataFrame
result_df = result.drop(result[rows_with_nan_or_blank].index)

result_df.shape

result_df.to_csv('features_column.csv', index=False)