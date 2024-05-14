import json
import requests


path = "../../source/pharmacodb/compounds_table.json"

def load_json(path):
        try:
        with open(path, encoding='utf-8') as f:
            this_json = json.load(f)
            return this_json
    except json.JSONDecodeError as e:
        print("Error decoding JSON: {}".format(e))


def fetch_pubchem_property(pubchem_id, compound_data):
    """ n ~ 54483 compounds"""
    response = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/property/MolecularFormula,MolecularWeight,,Fingerprint2D,InChIKey,InChI,CanonicalSMILES,IsomericSMILES/JSON")

    if response.status_code == 200:
        print(f"Successfully fetched {pubchem_id}")
        compound_dat = response.json()
        pubchem_properties = {"pubchem_properties": compound_dat['PropertyTable']['Properties'][0]}
        compound_data['annotation'].update(pubchem_properties)
        return compound_data
    else:
        print(f"Error: Unable to fetch data for  {pubchem_id}. Status code: {response.status_code}")
        return None


def update_compounds_pubchem_property(path)
    """compounds property annotations in pubchem"""
    dat = load_json(path)

    for compound in dat['data']['compounds']:
        if compound['annotation']['pubchem']:
            fetch_pubchem_property(compound['annotation']['pubchem'], compound)

    with open("../../source/pharmacodb/compounds_table_pubchem_updates.json", "w") as output_file:
        output_file.write(json.dumps(dat, indent=4))

update_compounds_pubchem_property(path)
