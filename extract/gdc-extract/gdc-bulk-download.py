#!/usr/bin/env python

import os
import sys
import json
import argparse

import requests
from subprocess import call
from pprint import pformat, pprint

URL_BASE="https://gdc-api.nci.nih.gov/v0/"


WORKFLOW_MAP = {
    "fpkm" : "HTSeq - FPKM",
    "cna" : "DNAcopy"
}

TYPE_MAP = {
    "clinical" : "clinical_supplement",
    "biospecimen" : "biospecimen_supplement",
    "methylation" : "methylation_beta_value"
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data_type")
    parser.add_argument("project_id")
    parser.add_argument("output_name")
    args = parser.parse_args()

    if args.data_type in WORKFLOW_MAP:
        query = {
            "op":"and",
            "content" : [
                {
                    "op":"=",
                    "content":{
                        "field":"analysis.workflow_type",
                        "value":[
                            WORKFLOW_MAP[args.data_type]
                        ]
                    }
                },
                {
                    "op":"=",
                    "content":{
                        "field":"cases.project.project_id",
                        "value":[
                            args.project_id
                        ]
                    }
                }
            ]
        }
    else:
        query = {
            "op":"and",
            "content" : [
                {
                    "op":"=",
                    "content":{
                        "field":"type",
                        "value":[
                            TYPE_MAP[args.data_type]
                        ]
                    }
                },
                {
                    "op":"=",
                    "content":{
                        "field":"cases.project.project_id",
                        "value":[
                            args.project_id
                        ]
                    }
                }
            ]
        }




    id_map = {}
    params = {}
    params['filters'] = json.dumps(query)
    params['expand'] = "cases.samples,cases.project"
    while 'size' not in params or data['pagination']['page'] < data['pagination']['pages']:
        params['size'] = 1000
        req = requests.get(URL_BASE + "files", params=params)
        data = req.json()['data']
        #print json.dumps(data, indent=4)
        for i in data['hits']:
            for case in i["cases"]:
                if "samples" in case:
                    for j in case["samples"]:
                        id_map[i['id']] = {"sample" : j['sample_id'], "project" : case["project"]["project_id"]}
                else:
                    id_map[i['id']] = {"project" : case["project"]["project_id"]}
        params['from'] = data['pagination']['from'] + data['pagination']['count']

    with open(args.output_name + ".map", "w") as handle:
        handle.write("file_id\tproject_id\tsample_id\n")
        for k, v in id_map.items():
            handle.write("%s\t%s\t%s\n" % (k, v["project"], v.get("sample", "")))

    print('downloading')
    headers = {'Content-type': 'application/json'}

    def chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i+n]

    keychunks = chunks(id_map.keys(), 100)
    paths = []
    for index, ids in enumerate(keychunks):
        r = requests.post(URL_BASE + 'data', data=json.dumps({"ids" : ids}), headers=headers, stream=True)
        path = args.output_name + str(index)
        paths.append(path)
        with open(path, 'wb') as f:
            for chunk in r.iter_content(1024):
                f.write(chunk)

    untar = "archive"

    call(["mkdir", "archive"])
    call(["mkdir", "manifest"])
    for index, path in enumerate(paths):
        call(["tar", "xzvf", path, "-C", "archive/"])
        call(["mv", "archive/MANIFEST.txt", "manifest/" + str(index)])

    call("cat manifest/* > archive/MANIFEST.txt", shell=True)

    tar = "cd " + untar + " && " + "tar czvf ../" + args.output_name + " . && cd .."
    print tar
    call(tar, shell=True)

    # call(["tar", "czvf", "../" + args.output_name, "."], cwd=untar)

    call(["rm", "-rf", untar])
    call(["rm", "-rf", "manifest"])
