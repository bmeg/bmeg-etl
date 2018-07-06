#!/usr/bin/env python

import sys
import pandas
import math
import json
from bmeg import phenotype_pb2, conversions, clinical_pb2
from google.protobuf import json_format

def proto_list_append(message, a):
    v = message.values.add()
    v.string_value = a

def gdsc_ic50_row(row, compound_table, sample_table, emit):
    sample_name, source = sample_table[row["COSMIC_ID"]]
    compound_name = compound_table[ int(row["DRUG_ID"]) ]

    gid = "gdsc-response:%s/%s" % (sample_name, compound_name)

    response = phenotype_pb2.ResponseCurve()
    # response.gid = gid
    response.source = source
    response.responseType = phenotype_pb2.ResponseCurve.ACTIVITY
    # response.compound = compound_name
    response.sample = sample_name

    compound = response.compounds.add()
    compound.ratio = 1.0
    compound.compound = compound_name

    s = response.summary.add()
    s.type = phenotype_pb2.ResponseSummary.IC50
    s.value = row['LN_IC50']
    s.unit = "uM"

    s = response.summary.add()
    s.type = phenotype_pb2.ResponseSummary.AUC
    s.value = row['AUC']
    s.unit = "uM"

    s = response.summary.add()
    s.type = phenotype_pb2.ResponseSummary.RMSE
    s.value = row['RMSE']
    s.unit = "uM"

    dose = row['CONC']
    """
    dr = response.values.add()
    dr.dose = dose
    dr.response = row['raw_max']
    for i in range(2,10):
        dose = dose / row['FOLD_DILUTION']
        dr = response.values.add()
        dr.dose = dose
        dr.response = row['raw%d' % (i) ]

    for i in range(1,49):
        v = row['control%d' % i]
        try:
            if not math.isnan(v):
                response.controls.append(v)
        except TypeError:
            pass
    for i in range(1,33):
        v = row['blank%d' % (i)]
        try:
            if not math.isnan(v):
                response.blanks.append(v)
        except TypeError:
            pass
    """
    emit(response)

def gdsc_cell_info(row, emit):
    sample = clinical_pb2.Biosample()
    sample.id = row["Sample Name"] # "gdsc:%s" % row["Sample Name"]
    sample.dataset_id = "gdsc"

    dis = row['GDSC\nTissue\ndescriptor 2']
    if not isinstance(dis, float):
        sample.disease.term = dis.lower().replace('_', ' ')

    sample.attributes['sampleType'] = "cellline"
    label = row['Cancer Type\n(matching TCGA label)']
    if not isinstance(label, float) and len(label):
        sample.attributes['source'] = label
    emit(sample)


class Emiter:
    def __init__(self, path):
        self.path = path
        self.handles = {}

    def emit(self, message):
        msg = json_format.MessageToDict(message)
        if message.DESCRIPTOR.full_name not in self.handles:
            self.handles[message.DESCRIPTOR.full_name] = open("%s.%s.json" % (self.path, message.DESCRIPTOR.full_name), "w")
        self.handles[message.DESCRIPTOR.full_name].write(json.dumps(msg) + "\n")

    def close(self):
        for v in self.handles.values():
            v.close()

conv_file = sys.argv[1]
cell_info_file = sys.argv[2]
compound_info_file = sys.argv[3]
raw_file = sys.argv[4]
fitted_file = sys.argv[5]
pubchem_file = sys.argv[6]


pubchem_table = {}
with open(pubchem_file) as handle:
    for line in handle:
        row = line.rstrip().split("\t")
        if row[1] == "none":
            pubchem_table[row[0]] = str(row[0]) # "compound:%s" % row[0]
        else:
            pubchem_table[row[0]] = str(row[1]) # "pubchem:%s" % row[1]

compound_table = {}
"""
comp_info = pandas.read_excel(conv_file, sheetname=1, index_col=0)
for row in comp_info.iterrows():
    compound_table[int(row[0])] = row[1]['GDSC name']
"""
comp_info = pandas.read_excel(compound_info_file, index_col=0)
for row in comp_info.iterrows():
    compound_table[int(row[0])] = conversions.pubchem(str(row[1]['Drug Name']))
    # compound_table[int(row[0])] = pubchem_table[str(row[1]['Drug Name'])]


e = Emiter("gdsc.scan")


cl_info = pandas.read_excel(conv_file, index_col=0)
sample_table = {}
for row in cl_info.iterrows():
    sample_table[row[0]] = (str(row[1]['CCLE name']), 'ccle') # "ccle:%s" % (row[1]['CCLE name'])

cl_info = pandas.read_excel(cell_info_file, index_col=1)
for row in cl_info.iterrows():
    if row[0] not in sample_table:
        sample_table[row[0]] = (str(row[1]['Sample Name']), 'gdsc') # "gdsc:%s" % (row[1]['Sample Name'])
        gdsc_cell_info(row[1], e.emit)

raw = pandas.read_csv(raw_file)
fitted = pandas.read_excel(fitted_file)

merge = pandas.merge(raw, fitted, on=["COSMIC_ID", "DRUG_ID"])

for r in merge.iterrows():
    cosmic_id = int(r[1]["COSMIC_ID"])
    if cosmic_id in cl_info.index:
        gdsc_ic50_row( r[1], compound_table, sample_table, e.emit )
