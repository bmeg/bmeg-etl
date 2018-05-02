#!/usr/bin/env python

import sys
import json
import requests

url = "https://api.gdc.cancer.gov/cases/"
headers = {'Content-Type': "application/json"}

base_query = {
    "expand": "samples.portions.analytes.aliquots",
    "fields": "samples.sample_id",
    # "from": 32500,
    # "from": 32750,
    "size": 250
}

output = open("tcga-barcode-gdc-uuid.tsv", 'a')
output.write("sample_id\taliquot_id\taliquot_barcode\n")

while True:
    try:
        r = requests.post(url, headers=headers, data=json.dumps(base_query))
        res = r.json()
        #print json.dumps(res, indent=4)
        for hit in res['data']['hits']:
            for sample in hit.get('samples', []):
                for portion in sample.get("portions", []):
                    for analyte in portion["analytes"]:
                        for aliquot in analyte["aliquots"]:
                            output.write("%s\t%s\t%s\n" % (sample["sample_id"], aliquot["aliquot_id"], aliquot["submitter_id"]))

        if res['data']["pagination"]["page"] < res['data']["pagination"]["pages"]:
            query["from"] = res['data']["pagination"]["from"] + res['data']["pagination"]["count"]
            print query["from"]

    except Exception:
        print(sys.exc_info()[0])

aliquot_query = {
    "expand": "samples",
    "fields": "samples.sample_id",
    "filters": {
        'op': 'in',
        'content': {
            'field': 'samples.portions.analytes.aliquots.submitter_id',
            'value': [sys.argv[1]]
        }
    }
}

# r = requests.post(url, headers=headers, data=json.dumps(aliquot_query))
# res = r.json()
# print json.dumps(res, indent=4)
# for hit in res['data']['hits']:
#     for sample in hit.get('samples', []):
#         print sample['sample_id']

