import json
import requests


def scrape_cases(page_ref=0):
    payload = {
        "format": "JSON",
        "expand": "annotations,demographic,diagnoses,diagnoses.treatments,exposures,family_histories,project,,samples,samples.portions,samples.portions.analytes,samples.portions.analytes.aliquots,tissue_source_site",
        "from": page_ref,
        "size": "1000"
    }

    response = requests.post("https://api.gdc.cancer.gov/cases", json=payload)
    response.raise_for_status()
    response_json = response.json()["data"]
    for hit in response_json["hits"]:
        print(json.dumps(hit))

    if "pagination" in response_json:
        if response_json["pagination"]["page"] < response_json["pagination"]["pages"]:
            page_ref = response_json["pagination"]["from"] + response_json["pagination"]["size"]
            scrape_cases(page_ref)
    return


if __name__ == "__main__":
    scrape_cases()
