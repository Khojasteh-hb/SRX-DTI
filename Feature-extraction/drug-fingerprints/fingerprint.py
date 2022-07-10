# Importing libraries
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys

def apply_getfinger(df_smile):
    finger_list = []
    for i in range(len(df_smile)):
        row = df_smile.iloc[i]
        mol = Chem.MolFromSmiles(row)
        fp = MACCSkeys.GenMACCSKeys(mol)
        finger_list.append(fp)
    return pd.Series(finger_list)

# Read compound ID of drugs
drug_smile = pd.read_csv('D:\my-projects\D-T prediction\Feature-extraction\drug-fingerprints\e_smile.csv')

# Remove unnamed columns
drug_smile.rename({"Unnamed: 0": "redundant"}, axis="columns", inplace=True)
drug_smile.drop(["redundant"], axis=1, inplace=True)

drug_finger = apply_getfinger(drug_smile['smile'])

drug_finger = pd.DataFrame(drug_finger, columns=['fingerprint_list'])

drug_finger['fingerprint'] = drug_finger.fingerprint_list.apply(lambda x: ', '.join([str(i) for i in x]))


# split column into multiple columns by delimiter
df_finger = pd.DataFrame(drug_finger['fingerprint'].str.split(',', expand=True))
#print(df_finger.head())

df_finger.to_csv('D:\my-projects\D-T prediction\Feature-extraction\drug-fingerprints\e_fingerprint.csv', index=True)







