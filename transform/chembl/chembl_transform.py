def build_citation(row):
    if "pubmed_id" in row.keys() and row["pubmed_id"]: 
        if "doi" in row.keys() and row["doi"]:
            cited_artifact = {
                "webLocation": [
                    {
                        "id": row["doi"],
                        "classifier": [
                            {
                                "coding": [
                                    {
                                        "system": "http://hl7.org/fhir/artifact-url-classifier",
                                        "code": "doi-based",
                                        "display": "DOI Based"
                                    }
                                ]
                            }
                        ],
                        "url": "https://doi.org/" + row["doi"]
                    }
                ]
            }
        else: 
            cited_artifact = None

        citation = {
            "resourceType": "Citation",
            "id": row["pubmed_id"],
            "identifier": [
                {
                    "system": "https://pubmed.ncbi.nlm.nih.gov/",
                    "value": row["pubmed_id"]
                }
            ],
            "status": "active",
            "citedArtifact": cited_artifact
        }
        return citation
    else: 
        pass


def build_substance_definition(row):

    if "synonyms" in row:
        synonymns = list(set(row["synonyms"]))
        synonyms_list = [{"name": s} for s in synonymns]
    else: 
        synonyms_list = [{"name": "synonym_name_place_holder"}]
        

    if "pref_name" in row and row["pref_name"]: 
        name = row["pref_name"]
    else: 
        name = "place_holder"

    if "pubmed_id" in row.keys() and row["pubmed_id"]:
        information_source = [{"reference": "Citation/" + row["pubmed_id"]}]
    else: 
        information_source = None

    if "morgan_fingerprint_2" in row: 
        classification = [{
            "coding": [
                            {
                                "system": "https://www.rdkit.org/",
                                "code": "Morgan Fingerprint Bit Vector",
                                "display": str(row["morgan_fingerprint_2"]) # use json.loads to revert back to list
                            }
                        ]
        }]
    else: 
        classification = None

    representation = []
    if "canonical_smiles" in row: 
        canonical_smiles_rep = {
                    "representation": row['canonical_smiles'],
                    "format": {
                        "coding": [
                            {
                                "system": "http://hl7.org/fhir/substance-representation-format",
                                "code": "SMILES",
                                "display": "SMILES"
                            }
                        ]
                    }
                }
        representation.append(canonical_smiles_rep)

    if "standard_inchi_key" in row:
        standard_inchi_key_rep = {
                    "representation": row['standard_inchi_key'],
                    "format": {
                        "coding": [
                            {
                                "system": "http://hl7.org/fhir/substance-representation-format",
                                "code": "InChI",
                                "display": "InChI"
                            }
                        ]
                    }
                }
        representation.append(standard_inchi_key_rep)

    if representation:
        structure = {
            "representation": representation
            }
    else: 
        structure = None


    substance_def = {
        "resourceType": "SubstanceDefinition",
        "id": row["chembl_id"],
        "name": [{"name": name, "synonym": synonyms_list }],
        "identifier": [
            {
                "system": "https://www.ebi.ac.uk/chembl/",
                "value": row["chembl_id"]
            }
        ],
        "informationSource": information_source,
        "classification": classification, 
        "structure": structure
    }
    # print(substance_def)
    return substance_def


def build_substance(row):
    substance = {
        "resourceType": "Substance",
        "id": row['chembl_id'],
        "instance": True,
        "code": {
            "reference": {
                "reference": "SubstanceDefinition/" + row['chembl_id']
            }
        }
    }
    return substance


def build_research_study(row):
    if "max_pahse" in row.keys() and row["max_pahse"]:
        phase = {
            "coding": [
                {
                    "system": "https://www.ebi.ac.uk/chembl/",
                    "code": "phase" + "-" + row["max_pahse"],
                    "display": "phase" + "-" + row["max_pahse"]
                }
            ]
        }
    else:
        phase = None

    research_study = {
        "resourceType": "ResearchStudy",
        "status": "unknown",
        "phase": phase,
        "focus": [
            {
                "reference": {
                    "reference": "Substance/" + row["chembl_id"]
                }
            }
        ]
    }
    return research_study
