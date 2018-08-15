#!/usr/bin/env python
"""
Bulk download legacy (hg19) data from the GDC

Edited on 2018-07-24 by Adam Struck
"""

import argparse
import json
import os
import re
import requests

from subprocess import call
from uuid import uuid4


URL_BASE = "https://api.gdc.cancer.gov/legacy/"

TYPE_MAP = {
    "expression": "Gene expression quantification",
    "cna": "Copy number segmentation",
    "methylation": "Methylation beta value",
    "clinical": "Clinical Supplement",
    "biospecimen": "Biospecimen Supplement"
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data_type")
    parser.add_argument("project_id")
    parser.add_argument("output_name")
    args = parser.parse_args()
    args.output_name = os.path.basename(args.output_name)

    query = {
        "op": "and",
        "content": [
            {
                "op": "=",
                "content": {
                    "field": "data_type",
                    "value": [
                        TYPE_MAP[args.data_type]
                    ]
                }
            },
            # Need these filters for cnv data
            # exlcude hg18 probes
            # exclude files with germline cnvs
            {
                "op": "exclude",
                "content": {
                    "field": "tags",
                    "value": [
                        "hg18",
                        "allcnv"
                    ]
                }
            },
            {
                "op": "=",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [
                        args.project_id
                    ]
                }
            }
        ]
    }

    data = {}
    id_map = {}
    params = {}
    params['filters'] = json.dumps(query)
    params['expand'] = "cases.samples,cases.project"

    print('creating file manifest...')

    while 'size' not in params or \
          data['pagination']['page'] < data['pagination']['pages']:
        params['size'] = 1000
        req = requests.get(URL_BASE + "files", params=params)
        data = req.json()['data']
        # print(json.dumps(data, indent=4))
        # sys.exit()
        for i in data['hits']:
            for case in i["cases"]:
                if "samples" in case:
                    for j in case["samples"]:
                        id_map[i['id']] = {
                            "sample": j['sample_id'],
                            "project": case["project"]["project_id"]
                        }
                else:
                    id_map[i['id']] = {
                        "project": case["project"]["project_id"]
                    }
        params['from'] = data['pagination']['from'] + \
            data['pagination']['count']

    with open(args.output_name + ".map", "w") as handle:
        handle.write("file_id\tproject_id\tsample_id\n")
        for k, v in id_map.items():
            handle.write("%s\t%s\t%s\n" %
                         (k, v["project"], v.get("sample", "")))

    def chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]

    headers = {'Content-type': 'application/json'}
    keychunks = chunks(list(id_map.keys()), 100)
    paths = []
    print('downloading files...')
    for index, ids in enumerate(keychunks):
        r = requests.post(URL_BASE + 'data',
                          data=json.dumps({"ids": ids}),
                          headers=headers,
                          stream=True)
        path = args.output_name + "-" + str(index)
        paths.append(path)
        with open(path, 'wb') as f:
            for chunk in r.iter_content(1024):
                f.write(chunk)

    print('creating archive...')

    random = str(uuid4())
    archive = "archive" + "-" + random
    manifest = "manifest" + "-" + random
    call(["mkdir", archive])
    call(["mkdir", manifest])
    for index, path in enumerate(paths):
        call(["tar", "xzvf", path, "-C", archive + "/"])
        call(["rm", "-f", path])
        call(["mv", archive + "/MANIFEST.txt", manifest + "/" + str(index)])

    output_name = re.sub("\.tar\.gz$", "", args.output_name) + ".tar.gz"
    call("cat" + manifest + "/* > " + archive + "/MANIFEST.txt", shell=True)
    tar = "cd " + archive + " && " + "tar czvf ../" + output_name + \
        " . && cd .."
    call(tar, shell=True)
    call(["rm", "-rf", archive])
    call(["rm", "-rf", manifest])
