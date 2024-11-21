#!/usr/bin/env python


def build_substance_definition(row):
    substance_def = {
        "resourceType": "SubstanceDefinition",
        "id": row['drug_id'],
        "identifier": [
            {
                "system": "https://pubchem.ncbi.nlm.nih.gov/",
                "value": row["pubchem_id"]
            }, 
            {
                "system": "https://pnnl-compbio.github.io/coderdata/",
                "value": row["improve_drug_id"]
            }
        ],
        "structure": {
            "representation": [
                {
                    "representation": row['canSMILES'],
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
                    "representation": row['InChIKey'],
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
        "id": row['drug_id'],
        "instance": True,
        "code": {
            "reference": {
                "reference": "SubstanceDefinition/" + row['drug_id']
            }
        }
    }
    return substance



def build_response(row):
    drug_response = {
        "aac": float(row['aac']),
        "auc": float(row['auc']),
        "dss1": float(row['dss']),
        "ec50": float(row['ec50']),
        "einf": float(row['einf']),
        "hs": float(row['hs']),
        "ic50": float(row['ic50']),
        "id": row['id'],
        "projectId": row['study'],
        "resourceType": "drug_response",
        "submitterId": row["improve_drug_id"] + "." + str(row["improve_sample_id"]),
        "substances": [
            {
                "code": {
                    "reference": {
                        "reference": "SubstanceDefinition/" +  row['drug_id'],
                        "resourceType": "Reference"
                    },
                    "resourceType": "CodeableReference"
                },
                "instance": True,
                "resourceType": "Substance"
            }
        ]
    }
    return drug_response


def build_patient(row):
    # Specimen has - sample.type info ex. cancer / normal
    substance = {
        "resourceType": "Substance",
        "id": row["drug_id"],
        "instance": True,
        "code": {
            "reference": {
                "reference": "SubstanceDefinition/" + row['drug_id']
            }
        }
    }

    patient = {"resourceType": "Patient",
               "id": row["sample_id"],
               "identifier": [{"system": "https://pnnl-compbio.github.io/coderdata/",
                               "value": row['id']}],
               "substances": [substance]}
    return patient

