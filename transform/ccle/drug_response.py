""" match drug response data to aliquot """

from types import SimpleNamespace as SN

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import DrugResponse, Aliquot
from bmeg.edge import ResponseIn, ResponseTo
from bmeg.enrichers.drug_enricher import compound_factory
from bmeg.ccle import build_ccle2depmap_conversion_table, missing_ccle_cellline_factory


NAMES = {
    'CCLE Cell Line Name': 'sample_id',
    # 'Primary Cell Line Name': 'primary_cell_line_name',
    'Compound': 'compound_id',
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


def transform(biosample_path='outputs/ccle/Biosample.Vertex.json.gz',
              drug_response_path='source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv',
              emitter_prefix='drug_response',
              emitter_directory="ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    # lookup table to convert to DepMap_IDs
    samples = build_ccle2depmap_conversion_table(biosample_path)

    # input and map
    input_stream = bmeg.ioutils.read_csv(drug_response_path)
    floats = ['amax', 'act_area', 'ec50', 'ic50', 'num_data']
    compound_gids = []
    # read the drug response csv
    missing_cell_lines = []
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

        # if no match, we will need to create project->individual->biosample->aliquot
        if drug_response.sample_id in samples:
            drug_response.sample_id = samples[drug_response.sample_id]
        elif drug_response.sample_id.split("_")[0] in samples:
            drug_response.sample_id = samples[drug_response.sample_id.split("_")[0]]
        else:
            if drug_response.sample_id not in missing_cell_lines:
                missing_cell_lines.append(drug_response.sample_id)

        # create drug_response vertex
        drug_resp = DrugResponse(**drug_response.__dict__)
        emitter.emit_vertex(drug_resp)
        #  and edge to aliquot
        emitter.emit_edge(
            ResponseIn(),
            drug_resp.gid(),
            Aliquot.make_gid('CCLE:{}'.format(drug_response.sample_id)),
        )
        # create compound
        compound = compound_factory(name=drug_response.compound_id)
        if compound.gid() not in compound_gids:
            emitter.emit_vertex(compound)
            compound_gids.append(compound.gid())
        # and an edge to it
        emitter.emit_edge(
            ResponseTo(),
            drug_resp.gid(),
            compound.gid(),
        )

    # generate project, individual, biosample, aliquot for missing cell lines
    missing_ccle_cellline_factory(emitter=emitter,
                                  source="CCLE",
                                  missing_ids=missing_cell_lines)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
