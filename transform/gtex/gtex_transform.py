

def build_patient(row):
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


