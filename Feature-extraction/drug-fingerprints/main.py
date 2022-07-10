# Importing libraries
import pandas as pd
from Pubdrug import getPubchem
from Pubdrug import getSmiles


def apply_getPubchem(df_ID):
    cid_list = []
    for i in range(len(df_ID)):
        row = df_ID.iloc[i]
        cid_list.append(getPubchem(row))
    return pd.Series(cid_list)


# Read compound ID of drugs
df_drug_ID = pd.read_csv('\E\e_drug_names.csv')
drug_ID = df_drug_ID['compound_ID']

cid_list = apply_getPubchem(drug_ID)
#print(cid_list.head())

drug_smile = getSmiles(cid_list)

drug_smile = pd.DataFrame(drug_smile, columns=['smile'])



df_smile = pd.concat([drug_ID, drug_smile], axis=1, sort=False)


df_smile.to_csv('\E\e_smile.csv', index=True)





