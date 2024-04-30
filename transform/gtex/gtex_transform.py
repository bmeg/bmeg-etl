

def build_patient(row):
    patient = {
        "resourceType": "Patient",
        "id": row["patient_id"],
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
    print(row)
    sample = {
        "resourceType": "Specimen",
        "id": row["sample_id"],
        "identifier": [
            {
                "system": "https://gtexportal.org/sample",
                "value": row["SAMPID"]
            }
        ],
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
            "reference": "Patient" + "/" + row["sample_patient_id"]
        }
    }
    return sample 


