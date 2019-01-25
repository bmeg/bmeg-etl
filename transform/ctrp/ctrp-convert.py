#!/usr/bin/env python

import os
import sys
import pandas
from bmeg.vertex import DrugResponse, Aliquot
from bmeg.edge import ResponseIn, ResponseTo

from bmeg.enrichers.drug_enricher import compound_factory

from bmeg.emitter import JSONEmitter

data_dir = sys.argv[1]

base_dir = os.path.dirname(__file__)

ccl_table_path = os.path.join(base_dir, "ctrp-cellline.table")
ccl_table = {}
with open(ccl_table_path) as handle:
    for line in handle:
        row = line.rstrip().split("\t")
        ccl_table[row[0]] = row[1]

metadrugPath = os.path.join(data_dir, "v20.meta.per_compound.txt")
metacelllinePath = os.path.join(data_dir, "v20.meta.per_cell_line.txt")
responsePath = os.path.join(data_dir, "v20.data.curves_post_qc.txt")
metaexperimentPath = os.path.join(data_dir, "v20.meta.per_experiment.txt")

curvePath = os.path.join(data_dir, "v20.data.per_cpd_post_qc.txt")

compound_df = pandas.read_table(metadrugPath)
compound_df = compound_df.set_index("master_cpd_id")

ccl_df = pandas.read_table(metacelllinePath)
ccl_df = ccl_df.set_index("master_ccl_id")

metaexperiment_df = pandas.read_table(metaexperimentPath)
metaexperiment_df = metaexperiment_df.drop_duplicates("experiment_id")
metaexperiment_df = metaexperiment_df.set_index("experiment_id")

response_df = pandas.read_table(responsePath)
curve_df = pandas.read_table(curvePath)
curve_df = curve_df.set_index(["experiment_id", "master_cpd_id"])

emitter = JSONEmitter(directory="ctrp")

for i, row in response_df.iterrows():
    exp_id = row['experiment_id']
    cpd_id = row['master_cpd_id']
    ccl_id = metaexperiment_df.loc[exp_id]['master_ccl_id']
    ccl_name = ccl_df.loc[ccl_id]['ccl_name']
    if ccl_name in ccl_table:
        ccl_name = ccl_table[ccl_name]
    cpd_name = compound_df.loc[cpd_id]['cpd_name']
    auc = row['area_under_curve']
    ec50 = row['apparent_ec50_umol']

    # curve_sub = curve_df.loc[ (curve_df["experiment_id"] == exp_id) & (curve_df["master_cpd_id"] == cpd_id) ]
    curve_sub = curve_df.loc[(exp_id, cpd_id)]
    conc = curve_sub["cpd_conc_umol"]
    resp = curve_sub["cpd_avg_pv"]

    dr = DrugResponse(sample_id=ccl_name, compound_id=cpd_name, source="CTRP",
                      act_area=auc, ec50=ec50, doses_um=list(conc),
                      activity_data_median=list(resp))
    compound = compound_factory(name=cpd_name)
    emitter.emit_vertex(dr)
    emitter.emit_edge(
        ResponseIn(),
        dr.gid(),
        Aliquot.make_gid(ccl_name)
    )
    emitter.emit_edge(
        ResponseTo(),
        dr.gid(),
        compound.gid()
    )


emitter.close()
