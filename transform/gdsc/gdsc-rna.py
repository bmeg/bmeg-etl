#!/usr/bin/env python

import sys
import json
import pandas
from bmeg import matrix_pb2
from google.protobuf import json_format
import numpy as np

# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/sanger1018_brainarray_ensemblgene_rma.txt.gz

def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    return json.dumps(msg)


"""
cl_info = pandas.read_excel(sys.argv[2], index_col=0)
ccle_sample_table = {}
gdsc_sample_table = {}
for row in cl_info.iterrows():
    ccle_sample_table[str(row[0])] = "biosample:CCLE:%s" % (row[1]['CCLE name'])
    gdsc_sample_table[str(row[0])] = "biosample:GDSC:%s" % (row[1]['GDSC1000 name'])
"""

cl_info = pandas.read_excel(sys.argv[2], index_col=0)
gdsc_sample_table = {}
for row in cl_info.iterrows():
    gdsc_sample_table[ "%0.f" % (row[1]["COSMIC identifier"]) ] = "gdsc:%s" % (row[0])

df = pandas.read_csv(sys.argv[1], sep="\t", index_col=0)
df = df[ np.invert(pandas.isnull(df.index)) ].drop('GENE_title', 1)

nc = []
for i in df.columns:
    name = i.split(".")[1]
    #if name in sample_table:
    #    nc.append(sample_table[name])
    #else:
    if name in gdsc_sample_table:
        nc.append(gdsc_sample_table[name])
    else:
        nc.append("biosample:cosmic:%s" % (i))

df.columns = pandas.Index(nc)
df = df.transpose()

for rid, row in df.iterrows():
    ge = matrix_pb2.GeneExpression()
    ge.gid = "GeneExpression:%s" % (rid)
    ge.biosample_id = rid
    for k, v in row.iteritems():
        ge.expressions[k] = v
    
    print message_to_json(ge)
    #handle.write("%s\n" % (message_to_json(ge)))