

def build_patient(row):
    """
    import orjson
    from fhir.resources.patient import Patient
    from fhir.resources.identifier import Identifier

    patient_identifier = Identifier(**{"value": row['SUBJID'], "system": "".join(["https://gtexportal.org/", "subject"])})
    patient = Patient(**{"id": "-".join([row['SUBJID'], "gtex"]), "identifier": [patient_identifier], "birthDate": row['fhir_birthDate'], "gender": row['SEX_string']})
    return orjson.loads(patient.json())
    """

    # TODO: replace and confirm BMEG id def
    patient = {
        "resourceType": "Patient",
        "id": row["SUBJID"] + "-gtex",
        "identifier": [
            {
                "system": "https://gtexportal.org/subject",
                "value": row["SUBJID"]
            }
        ],
        "gender": row["SEX_string"],
        "birthDate": row["fhir_birthDate"]
    }
    return patient



def build_specimen(row):
    """
    import orjson
    from fhir.resources.specimen import Specimen
    from fhir.resources.reference import Reference
    from fhir.resources.codeableconcept import CodeableConcept

    cc = CodeableConcept(**{"coding": [{"system": "http://snomed.info/sct", "display": row['SMTSD'], "code": str(row['sctid'])}]})
    subject_ref = Reference(**{"reference": "/".join(["Patient", "-".join([row['SUBJID'], "gtex"])])})
    specimen = Specimen(**{"id": "-".join([row['SAMPID'], "gtex"]), "type": cc, "subject": subject_ref})
    return orjson.loads(specimen.json())
    """
    # TODO: replace and confirm BMEG id def
    sample = {
        "resourceType": "Specimen",
        "id": row["SAMPID"] + "-gtex",
        "type": {
            "coding": [
                {
                    "system": "http://snomed.info/sct",
                    "code": str(row["sctid"]),
                    "display": row["SMTSD"]
                }
            ]
        },
        "subject": {
            "reference": "Patient" + "/" + row["SUBJID"] + "-gtex"
        }
    }
    return sample 


