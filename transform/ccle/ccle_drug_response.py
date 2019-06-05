""" match drug response data to aliquot """

from types import SimpleNamespace as SN

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import (Aliquot, DrugResponse, Project,
                  Compound_Projects_Project,
                  DrugResponse_Aliquot_Aliquot,
                  DrugResponse_Compound_Compound)
from bmeg.enrichers.drug_enricher import compound_factory


NAMES = {
    'CCLE Cell Line Name': 'submitter_id',
    # 'Primary Cell Line Name': 'primary_cell_line_name',
    'Compound': 'submitter_compound_id',
    # 'Target': 'target',
    'Doses (uM)': 'doses_um',
    'Activity Data (median)': 'activity_data_median',
    'Activity SD': 'activity_sd',
    'Num Data': 'num_data',
    'FitType': 'fit_type',
    'EC50 (uM)': 'ec50',
    'IC50 (uM)': 'ic50',
    'Amax': 'amax',
    'ActArea': 'act_area',
}


def transform(cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              drug_response_path='source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv',
              emitter_prefix='drug_response',
              emitter_directory="ccle"):

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
    projects = bmeg.ioutils.read_lookup(project_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    # input and map
    input_stream = bmeg.ioutils.read_csv(drug_response_path)
    floats = ['amax', 'act_area', 'ec50', 'ic50', 'num_data']
    compound_gids = {}

    # read the drug response csv
    project_compounds = {}
    for line in input_stream:
        # map the names to snake case
        mline = {"source": "CCLE"}
        for k, v in NAMES.items():
            # line[v] = line.pop(k)
            # csv field to float array
            if ',' in line[k]:
                mline[v] = [float(s) for s in line[k].split(',')]
            # NA -> None
            elif line[k] == 'NA':
                mline[v] = None
            # string scalar to float
            elif line[k] and v in floats:
                mline[v] = float(line[k])
            else:
                mline[v] = line[k]
        drug_response = SN(**mline)

        i = drug_response.submitter_id
        if i in celllines:
            cellline_id = celllines[i]
        elif i.split("_")[0] in celllines:
            cellline_id = celllines[i.split("_")[0]]
        else:
            cellline_id = i
        drug_response.submitter_id = DrugResponse.make_gid(cellline_id)

        project_id = "CCLE_%s" % (projects.get(cellline_id, "Unknown"))
        proj = Project(submitter_id=Project.make_gid(project_id), project_id=project_id)
        if proj.gid() not in project_compounds:
            project_compounds[proj.gid()] = {}

        # create drug_response vertex
        drug_resp = DrugResponse(**drug_response.__dict__)
        emitter.emit_vertex(drug_resp)
        #  and edge to aliquot
        emitter.emit_edge(
            DrugResponse_Aliquot_Aliquot(
                from_gid=drug_resp.gid(),
                to_gid=Aliquot.make_gid("CCLE:%s:DrugResponse:%s" % (drug_response.submitter_id, drug_response.submitter_compound_id))
            ),
            emit_backref=True
        )
        # create compound
        compound = compound_factory(name=drug_response.submitter_compound_id)
        if compound.gid() not in compound_gids:
            emitter.emit_vertex(compound)
            compound_gids[compound.gid()] = None
        # and an edge to it
        emitter.emit_edge(
            DrugResponse_Compound_Compound(
                from_gid=drug_resp.gid(),
                to_gid=compound.gid()
            ),
            emit_backref=True
        )

        # create edge from compound to project
        if compound.gid() not in project_compounds[proj.gid()]:
            emitter.emit_edge(
                Compound_Projects_Project(
                    from_gid=compound.gid(),
                    to_gid=proj.gid()
                ),
                emit_backref=True
            )
            project_compounds[proj.gid()][compound.gid()] = True

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
