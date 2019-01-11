""" match drug response data to aliquot """

import glob
import ujson
from types import SimpleNamespace as SN

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import DrugResponse, Aliquot  # , Biosample, Individual, Project
from bmeg.edge import ResponseIn, ResponseTo  # , AliquotFor, BiosampleFor, InProject
from bmeg.enrichers.drug_enricher import compound_factory


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


def transform(biosample_path='outputs/ccle/Biosample.Vertex.json*',
              drug_response_path='source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv',
              emitter_prefix='drug_response',
              emitter_directory="ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

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
    floats = ['amax', 'act_area', 'ec50', 'ic50', 'num_data']
    compound_gids = []
    # read the drug response csv
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

        sample_id = drug_response.sample_id
        # create drug_response vertex
        pharmacalogical_profile = DrugResponse(**drug_response.__dict__)
        emitter.emit_vertex(pharmacalogical_profile)
        #  and edge to aliquot
        emitter.emit_edge(
            ResponseIn(),
            pharmacalogical_profile.gid(),
            Aliquot.make_gid(sample_id),
        )
        # create compound
        compound = compound_factory(name=drug_response.compound_id)
        if compound.gid() not in compound_gids:
            emitter.emit_vertex(compound)
            compound_gids.append(compound.gid())
        # and an edge to it
        emitter.emit_edge(
            ResponseTo(),
            pharmacalogical_profile.gid(),
            compound.gid(),
        )

    """
    # create any missing vertexes
    individual_gids = project_gids = []
    for ccle_id in missing_cell_lines:
        b = Biosample(ccle_id)
        emitter.emit_vertex(b)
        a = Aliquot(aliquot_id=ccle_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid(),
        )
        i = Individual(individual_id='CCLE:{}'.format(ccle_id))
        if i.gid() not in individual_gids:
            emitter.emit_vertex(i)
            individual_gids.append(i.gid())
        emitter.emit_edge(
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
            emitter.emit_vertex(p)
            project_gids.append(p.gid())
        emitter.emit_edge(
            InProject(),
            i.gid(),
            p.gid(),
        )
    """

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
