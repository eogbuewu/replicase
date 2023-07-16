# Install the ChEMBL web service package so that we can retrieve bioactivity data from the ChEMBL Database.
! pip install chembl_webresource_client

# Import necessary libraries
import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client

# Target search for coronavirus
target = new_client.target
target_query = target.search('coronavirus')
targets = pd.DataFrame.from_dict(target_query)
targets

#We will assign the 8,9 (Replicase polyprotein 1ab) to the *selected_target* variable 
selected_target1 = targets.target_chembl_id[8] #binding assays (CHEMBL5118)
selected_target2 = targets.target_chembl_id[9] #functional assays (CHEMBL4523582)
print("selected_target1 :", selected_target1)
print("selected_target2 :", selected_target2)

# Here, we will retrieve only bioactivity data for coronavirus replicase polyprotein 1ab that are reported as IC
# values in nM (nanomolar) unit.
activity = new_client.activity
res1 = activity.filter(target_chembl_id=selected_target1).filter(standard_type="IC50")
res2 = activity.filter(target_chembl_id=selected_target2).filter(standard_type="IC50")

#put them in a dataframe 
df1 = pd.DataFrame.from_dict(res1)
df2 = pd.DataFrame.from_dict(res2)

#Concatenate the dataframes vertically
data = pd.concat([df1, df2], axis=0)
# Reset the index of the resulting dataframe
data = data.reset_index(drop=True)
data.head(5)
data.shape

#Finally we will save the resulting bioactivity data to a CSV file bioactivity_data.csv.
data.to_csv('replicase_data_raw.csv', index=False)

# If any compounds has missing value for the standard_value column then drop it
mydata = data[data.standard_value.notna()]
mydata
mydata.shape
#handle compounds that have no smiles.
newdata = mydata[mydata.canonical_smiles.notna()]
newdata.shape
#check to see if there are duplicates
len(newdata.canonical_smiles.unique())
#Drop the duplicates
newdata = newdata.drop_duplicates(['canonical_smiles'])
newdata.shape
#Reset the index
newdata.reset_index(drop=True, inplace=True)

#calculate the pChEMBL value of the rows that are missing using the standard value of the IC50 
newdata['standard_value'] = pd.to_numeric(newdata['standard_value'], errors='coerce')
newdata['pchembl_value'] = pd.to_numeric(newdata['pchembl_value'], errors='coerce')
newdata.loc[newdata['pchembl_value'].isnull(), 'pchembl_value'] = np.log10(newdata.loc[newdata['pchembl_value'].isnull(), 'standard_value'].values)
newdata['pchembl_value'] = newdata['pchembl_value'].round(2)

bioactivity_class = []
for i in newdata.pchembl_value:
    if float(i) <= 4:
        bioactivity_class.append("no activity")
    elif float(i) > 4 and float(i) <= 5.99:
        bioactivity_class.append("low activity")
    elif float(i) >= 6 and float(i) <= 7.99:
        bioactivity_class.append("moderate activity")
    else:
        bioactivity_class.append("high activity")

len(bioactivity_class)

#slice your data
mynewdata= newdata[['molecule_chembl_id','canonical_smiles','pchembl_value']]
mynewdata.shape
#serialize the bioactivity list
bioactivity_class = pd.Series(bioactivity_class, name='bioactivity_class')

# Reset the index of 'mynewdata' and 'bioactivity_class' to ensure consistent indices
mynewdata.reset_index(drop=True, inplace=True)
bioactivity_class.reset_index(drop=True, inplace=True)

# Concatenate 'mynewdata' and 'bioactivity_class' horizontally (axis=1)
dataframe = pd.concat([mynewdata, bioactivity_class], axis=1)

# Verify the resulting DataFrame
print(dataframe)
print(dataframe)
#save the preprocessed data
dataframe.to_csv('Replicase_bioactivity_data_preprocessed.csv', index=False)

       
    
    
 

    

