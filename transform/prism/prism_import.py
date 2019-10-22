#!/usr/bin/env python

import os
import sys
import pandas

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import (Aliquot, DrugResponse, Project, Compound,
                  Compound_Projects_Project,
                  DrugResponse_Aliquot_Aliquot,
                  DrugResponse_Compounds_Compound)
from bmeg.enrichers.drug_enricher import compound_factory

base = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    datadir = os.path.abspath(sys.argv[1])
    #mfi = pandas.read_csv(os.path.join(datadir, "primary_MFI.csv"), index_col=0)

    compound_map = {}
    with open(os.path.join(base, "compound.table")) as handle:
        for line in handle:
            row = line.rstrip().split("\t")
            compound_map[row[0]] = row[1]

    transformed = pandas.read_csv(os.path.join(datadir, "primaryscreenreplicatecollapsedlogfoldchange.csv"), index_col=0)
    treatment = pandas.read_csv(os.path.join(datadir, "primaryscreenreplicatecollapsedtreatmentinfo.csv"), index_col=0)

    emitter = JSONEmitter(prefix="drug_response", directory="prism")

    emitted_compounds = {}
    project_compounds = {}
    for rowname, row in transformed.iterrows():
        for colname, val in row.items():

            brd_name = treatment.broad_id[colname]
            sub_name = treatment.name[colname]
            if brd_name in compound_map:
                drug_id = compound_map[drug_name]
            elif sub_name in compound_map:
                drug_id = compound_map[sub_name]
            else:
                print("Missing Compound: %s" % (sub_name))
                drug_id = "TODO:%s" % (brd_name)

            cellline_id = rowname
            # create drug response vertex
            dr = DrugResponse(id=DrugResponse.make_gid("PRISM", rowname, colname, drug_id),
                              aac=val,
                              submitter_id=colname,
                              #doses_um=[ treatment.dose[colname] ],
                              project_id="PRISM")
            emitter.emit_vertex(dr)

            emitter.emit_edge(
                DrugResponse_Aliquot_Aliquot(
                    from_gid=dr.gid(),
                    to_gid=Aliquot.make_gid("PRISM:%s" % (cellline_id))
                ),
                emit_backref=True
            )

            # create compound
            #compound = compound_factory(name=drug_name)
            compound = Compound(id=Compound.make_gid(drug_id),
                            id_source='TODO',
                            submitter_id=brd_name,
                            project_id=Project.make_gid('Reference'))

            #if compound.gid() not in emitted_compounds:
            #    emitter.emit_vertex(compound)
            #    emitted_compounds[compound.gid()] = True

            emitter.emit_edge(
                DrugResponse_Compounds_Compound(
                    from_gid=dr.gid(),
                    to_gid=compound.gid()
                ),
                emit_backref=True
            )

    emitter.close()
