def build_specimen(row):
    if "Cellosaurus.Accession.id" in row.keys() and row["Cellosaurus.Accession.id"]:
        subject = {
            "reference": "Patient" + "/" + row["Cellosaurus.Accession.id"]
        }
    else: 
        subject = None

    sample = {
        "resourceType": "Specimen",
        "id": row["sampleid"],
        "identifier": [
            {
                "system": "https://pharmacodb.ca/",
                "value": row["sampleid"]
            }
        ],
        "subject": subject
    }
    return sample 