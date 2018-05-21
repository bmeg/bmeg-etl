#!/usr/bin/env python

import sys
import yaml
import json

with open(sys.argv[1]) as handle:
    data = yaml.load(handle.read())

vertices = [{"gid" : "root", "label" : "Query"}]
edges = []
for ent in data:
    for i in ent.get('vertexes', []):
        o = {"gid" : i['label'], "label" : "Object", "data" : {'fields':{}}}
        for f in i.get('fields', []):
            o['data']['fields'][f] = "string"
        vertices.append(o)

    for i in ent.get("edges", []):
        o = {"from" : i['fromLabel'], "to" : i["toLabel"], "label" : "field",
            "data" : {
                "label" : i["label"],
                "name" : i["label"]
            }
        }
        edges.append(o)
graph = {
    "vertices" : vertices,
    "edges" : edges
}

with open("bmeg.schema.nodes", "w") as handle:
    for i in graph["vertices"]:
        handle.write(json.dumps(i) + "\n")

with open("bmeg.schema.edges", "w") as handle:
    for i in graph["edges"]:
        handle.write(json.dumps(i) + "\n")
