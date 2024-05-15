import json
import requests


def load_json(path):
    try:
        with open(path, encoding='utf-8') as f:
            this_json = json.load(f)
            return this_json
    except json.JSONDecodeError as e:
        print("Error decoding JSON: {}".format(e))


def make_request(api_url, retries=3):
    delay = 0.5
    for _ in range(retries):
        response = requests.get(api_url)

        if response.status_code == 200:
            return response.json()
        else:
            print(f"Received status code: {response.status_code}. Retrying...")
            delay *= 2 ** retries  # change delay
            time.sleep(delay + random.uniform(0, 1))  # add jitter
    raise Exception("Failed to fetch data after multiple retries")


def fetch_pubchem_property(pubchem_id, compound_data):
    """ n ~ 54483 compounds"""
    api_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/property/MolecularFormula,MolecularWeight,Fingerprint2D,InChIKey,InChI,CanonicalSMILES,IsomericSMILES/JSON"

    try:
        compound_dat = make_request(api_url, retries=3)
        print(f"Successfully fetched {pubchem_id}")
        pubchem_properties = {"pubchem_properties": compound_dat['PropertyTable']['Properties'][0]}
        compound_data['annotation'].update(pubchem_properties)
        return compound_data
    except Exception as e:
        print(f"Error: Unable to fetch data for {pubchem_id}: {e}")
        return None


def update_compounds_pubchem_property(path):
    """compounds property annotations in pubchem"""
    dat = load_json(path)

    for compound in dat['data']['compounds']:
        if compound['annotation']['pubchem']:
            fetch_pubchem_property(compound['annotation']['pubchem'], compound)

    # with open("../../source/pharmacodb/graphql_dump/compounds_table_pubchem_updates.json", "w") as output_file:
    #    output_file.write(json.dumps(dat, indent=4))

    # save as ndjson 
    with open("../../source/pharmacodb/graphql_dump/compounds_table_pubchem_updates.ndjson", 'w', encoding='utf8') as file:
        file.write('\n'.join(map(lambda e: json.dumps(e, ensure_ascii=False), dat['data']['compounds'])))


update_compounds_pubchem_property(path="../../source/pharmacodb/graphql_dump/compounds_table.json")


