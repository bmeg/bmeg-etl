def build_substance_definition(row):
    # "informationSource": Reference(Citation) for publications
    substance_def = {
        "resourceType": "SubstanceDefinition",
        "id": row["chembl_id"],
        "name": [{"name": row["pref_name"], "synonym": [{"name": "synonym_name_place_holder"}] }],
        "identifier": [
            {
                "system": "https://www.ebi.ac.uk/chembl/",
                "value": row["chembl_id"]
            }
        ],
        "structure": {
            "representation": [
                {
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
                },
                {
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
            ]
        }
    }
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
    # ResearchStudy.phase -> codeableConcept 
    # ResearchStudy.focus -> Substance
    return 
