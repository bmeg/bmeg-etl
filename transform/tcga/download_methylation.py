import argparse
import json
import requests
import sys

from subprocess import call
from uuid import uuid4


URL_BASE = "https://api.gdc.cancer.gov/legacy/"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--platform",
        choices=["Illumina Human Methylation 450",
                 "Illumina Human Methylation 27"],
        default="Illumina Human Methylation 450",
        help="download data from which methylation platform"
    )
    parser.add_argument(
        "--manifest-only",
        default=False,
        action="store_true",
        help="write file manifest, but do not download the files"
    )
    args = parser.parse_args()
    output_name = "source/tcga/methylation/" + args.platform.replace(" ", "")

    query = {
        "op": "and",
        "content": [
            {
                "op": "=",
                "content": {
                    "field": "data_type",
                    "value": [
                        "Methylation beta value"
                    ]
                }
            },
            {
                "op": "=",
                "content": {
                    "field": "cases.project.program.name",
                    "value": [
                        "TCGA"
                    ]
                }
            },
            {
                "op": "=",
                "content": {
                    "field": "platform",
                    "value": [
                        args.platform
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

    print('querying GDC API...')

    while 'size' not in params or \
          data['pagination']['page'] < data['pagination']['pages']:

        params['size'] = 1000
        req = requests.get(URL_BASE + "files", params=params)
        data = req.json()['data']
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

        print('processed page', data['pagination']['page'])

    print('creating file manifest...')
    output_manifest = output_name + ".map"
    with open(output_manifest, "w") as handle:
        handle.write("file_id\tproject_id\tsample_id\n")
        for k, v in id_map.items():
            handle.write("%s\t%s\t%s\n" %
                         (k, v["project"], v.get("sample", "")))

    if args.manifest_only:
        sys.exit(0)

    def chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]

    headers = {'Content-type': 'application/json'}
    keychunks = chunks(list(id_map.keys()), 100)
    paths = []
    print('downloading files...')
    headers = {'Content-type': 'application/json'}
    for index, ids in enumerate(keychunks):
        r = requests.post(URL_BASE + 'data',
                          data=json.dumps({"ids": ids}),
                          headers=headers,
                          stream=True)
        path = output_name + "-" + str(index)
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

    output_archive = output_name + ".tar.gz"
    call("cat" + manifest + "/* > " + archive + "/MANIFEST.txt", shell=True)
    tar = "cd " + archive + " && " + "tar czvf ../" + output_archive + \
        " . && cd .."
    call(tar, shell=True)
    call(["rm", "-rf", archive])
    call(["rm", "-rf", manifest])
