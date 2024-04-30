#!/usr/bin/env python


def build_substance_definition(row):
    substance_def = {
        "resourceType": "SubstanceDefinition",
        "id": str(int(row['pubchem_id'])),
        "identifier": [
            {
                "system": "https://pubchem.ncbi.nlm.nih.gov/",
                "value": str(int(row['pubchem_id']))
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
        "instance": True,
        "code": {
            "reference": {
                "reference": "SubstanceDefinition/" + str(int(row['pubchem_id']))
            }
        }
    }
    return substance



def build_response(row):
    drug_response = {
        "aac": row['aac'],
        "auc": row['auc'],
        "dss1": row['dss1'],
        "ec50": row['ec50'],
        "einf": row['einf'],
        "hs": row['hs'],
        "ic50": row['ic50'],
        "id": row['id'],
        "projectId": row['id'],
        "resourceType": "drug_response",
        "submitterId": row["improve_drug_id"] + ":" + str(row["improve_sample_id"]),
        "substances": [
            {
                "code": {
                    "reference": {
                        "reference": "SubstanceDefinition/" + str(int(row["pubchem_id"])),
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
        "instance": True,
        "code": {
            "reference": {
                "reference": "SubstanceDefinition/" + str(int(row["pubchem_id"]))
            }
        }
    }

    patient = {"resourceType": "Patient",
               "id": str(row["improve_sample_id"]) + ":" + row["improve_drug_id"],
               "identifier": [{"system": "https://pnnl-compbio.github.io/coderdata/",
                               "value": str(row['improve_sample_id'])}],
               "substances": [substance]}
    return patient

