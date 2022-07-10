# Importing libraries
import requests

# Get PubChem CID by compound ID
def getPubchem(compound_ID):
    """
    Get PubChem CID from the Compound IDs.

    Parameters
    ----------
    compound ID : string

    Returns
    -------
    PubChem CID : string
    """
    name = compound_ID
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"

    r = requests.get(url)
    r.raise_for_status()
    response = r.json()
    if "IdentifierList" in response:
        cid = response["IdentifierList"]["CID"][0]
    else:
        raise ValueError(f"Could not find matches for compound: {name}")
    return cid



# Get smile by PubChem CID
def getSmiles(cids):
    """
    Get the canonical SMILES string from the PubChem CIDs.

    Parameters
    ----------
    cids : list
        A list of PubChem CIDs.

    Returns
    -------
    list
        The canonical SMILES strings of the PubChem CIDs.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(map(str, cids))}/property/CanonicalSMILES/JSON"
    r = requests.get(url)
    r.raise_for_status()
    return [item["CanonicalSMILES"] for item in r.json()["PropertyTable"]["Properties"]]