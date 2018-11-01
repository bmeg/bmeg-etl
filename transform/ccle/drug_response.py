""" match drug response data to aliquot """

import logging
import glob
import ujson
from types import SimpleNamespace as SN

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import ParamacalogicalProfile, Aliquot, Biosample, Individual, Project
from bmeg.edge import ParamacalogicalProfileIn, ResponseTo, AliquotFor, BiosampleFor, InProject
from bmeg.enrichers.drug_enricher import compound_factory


NAMES = {
    'CCLE Cell Line Name': 'ccle_cell_line_name',
    'Primary Cell Line Name': 'primary_cell_line_name',
    'Compound': 'compound',
    'Target': 'target',
    'Doses (uM)': 'doses_um',
    'Activity Data (median)': 'activity_data_median',
    'Activity SD': 'activity_sd',
    'Num Data': 'num_data',
    'FitType': 'fit_type',
    'EC50 (uM)': 'ec50_um',
    'IC50 (uM)': 'ic50_um',
    'Amax': 'a_max',
    'ActArea': 'act_area',
}


def transform(biosample_path='outputs/ccle/Biosample.Vertex.json*',
              drug_response_path='source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv',
              emitter_prefix=None,
              emitter_directory="ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    # we create new [Aliquot, Biosample, Individual, Project] vertexes for items not in samples
    # A separate emitter, with its own prefix is used as to not overwrite existing vertexes and edges created by samples
    xtra_emitter = JSONEmitter(directory=emitter_directory, prefix='drug_response')

    # lookup table of Broad_ID
    samples = {}
    for path in glob.iglob(biosample_path):  # match .json or .json.gz
        input_stream = bmeg.ioutils.reader(path)
        for line in input_stream:
            biosample = SN(**ujson.loads(line))
            ccle_attributes = SN(**biosample.data['ccle_attributes'])
            # mangle the lookup keys
            samples[ccle_attributes.CCLE_Name] = ccle_attributes.Broad_ID
            samples[ccle_attributes.CCLE_Name.split('_')[0]] = ccle_attributes.Broad_ID
            samples[ccle_attributes.Aliases] = ccle_attributes.Broad_ID
        break
    assert len(samples.keys()), 'No biosamples found, does {} exist?'.format(biosample_path)

    # input and map
    input_stream = bmeg.ioutils.read_csv(drug_response_path)
    floats = ['a_max', 'act_area', 'ec50_um', 'ic50_um', 'num_data']
    missing_cell_lines = []
    compound_gids = []
    # read the drug response csv
    for line in input_stream:
        # map the names to snake case
        for k, v in NAMES.items():
            line[v] = line.pop(k)
            # csv field to float array
            if ',' in line[v]:
                line[v] = [float(s) for s in line[v].split(',')]
            # NA -> None
            if line[v] == 'NA':
                line[v] = None
            # string scalar to float
            if line[v] and v in floats:
                line[v] = float(line[v])
        drug_response = SN(**line)
        # match to existing biosample/aliquot
        broad_id = None
        # mangle the lookup keys
        keys = [
            drug_response.ccle_cell_line_name,
            drug_response.ccle_cell_line_name.split('_')[0],
            drug_response.primary_cell_line_name,
            drug_response.primary_cell_line_name.replace('-', ''),
        ]
        for k in keys:
            if k in samples:
                broad_id = samples[k]
        # if no match, we will need to create project->individual->biosample->aliquot
        sample_id = broad_id

        if not broad_id:
            sample_id = drug_response.ccle_cell_line_name
            if drug_response.ccle_cell_line_name not in missing_cell_lines:
                logging.debug('no-match', drug_response.ccle_cell_line_name, drug_response.primary_cell_line_name)
                missing_cell_lines.append(drug_response.ccle_cell_line_name)

        # create drug_response vertex
        pharmacalogical_profile = ParamacalogicalProfile(**drug_response.__dict__)
        emitter.emit_vertex(pharmacalogical_profile)
        #  and edge to aliquot
        emitter.emit_edge(
            ParamacalogicalProfileIn(),
            pharmacalogical_profile.gid(),
            Aliquot.make_gid(sample_id),
        )
        # create compound
        compound = compound_factory(name=drug_response.compound)
        if compound.gid() not in compound_gids:
            emitter.emit_vertex(compound)
            compound_gids.append(compound.gid())
        # and an edge to it
        emitter.emit_edge(
            ResponseTo(),
            pharmacalogical_profile.gid(),
            compound.gid(),
        )

    # create any missing vertexes
    individual_gids = project_gids = []
    for ccle_id in missing_cell_lines:
        b = Biosample(ccle_id)
        xtra_emitter.emit_vertex(b)
        a = Aliquot(aliquot_id=ccle_id)
        xtra_emitter.emit_vertex(a)
        xtra_emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid(),
        )
        i = Individual(individual_id='CCLE:{}'.format(ccle_id))
        if i.gid() not in individual_gids:
            xtra_emitter.emit_vertex(i)
            individual_gids.append(i.gid())
        xtra_emitter.emit_edge(
            BiosampleFor(),
            b.gid(),
            i.gid(),
        )
        # first see if we have wholesale name changes
        project_id = ccle_id
        # strip off prefix
        name_parts = project_id.split('_')
        name_start = 1
        if len(name_parts) == 1:
            name_start = 0
        project_id = '_'.join(name_parts[name_start:])
        # create project
        p = Project(project_id='CCLE:{}'.format(project_id))
        if p.gid() not in project_gids:
            xtra_emitter.emit_vertex(p)
            project_gids.append(p.gid())
        xtra_emitter.emit_edge(
            InProject(),
            i.gid(),
            p.gid(),
        )

    emitter.close()
    xtra_emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
