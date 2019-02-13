
import os
import contextlib
import pytest
import json
from transform.gdc.cases import transform
from bmeg.vertex import Sample, Aliquot, Case, Project, Compound
from bmeg.emitter import JSONEmitter
from bmeg.ioutils import reader


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/')


def validate(helpers, emitter_path_prefix, parameters):
    """ run xform and test results"""
    sample_file = '{}Sample.Vertex.json.gz'.format(emitter_path_prefix)
    aliquot_file = '{}Aliquot.Vertex.json.gz'.format(emitter_path_prefix)
    case_file = '{}Case.Vertex.json.gz'.format(emitter_path_prefix)
    compound_file = '{}Compound.Vertex.json.gz'.format(emitter_path_prefix)

    samplefor_file = '{}SampleFor.Edge.json.gz'.format(emitter_path_prefix)
    aliquotfor_file = '{}AliquotFor.Edge.json.gz'.format(emitter_path_prefix)
    inproject_file = '{}InProject.Edge.json.gz'.format(emitter_path_prefix)
    treatedwith_file = '{}TreatedWith.Edge.json.gz'.format(emitter_path_prefix)

    all_files = [sample_file, aliquot_file, case_file, samplefor_file,
                 aliquotfor_file, inproject_file, treatedwith_file, compound_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    emitter = JSONEmitter(emitter_path_prefix)
    transform(emitter, parameters)
    emitter.close()
    # test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Sample, sample_file)
    # test.Aliquot.Vertex.json
    helpers.assert_vertex_file_valid(Aliquot, aliquot_file)
    # test.Case.Vertex.json
    helpers.assert_vertex_file_valid(Case, case_file)
    # test.Compound.Vertex.json
    helpers.assert_vertex_file_valid(Compound, compound_file)
    # test.SampleFor.Edge.json
    helpers.assert_edge_file_valid(Sample, Case, samplefor_file)
    # test.AliquotFor.Edge.json
    helpers.assert_edge_file_valid(Aliquot, Sample, aliquotfor_file)
    # test.InProject.Edge.json
    helpers.assert_edge_file_valid(Case, Project, inproject_file)
    # test.TreatedWith.Edge.json
    helpers.assert_edge_file_valid(Case, Compound, treatedwith_file)

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Project']
    )

    # test Aliquot contents
    with reader(aliquot_file) as f:
        for line in f:
            # should be json
            aliquot = json.loads(line)
            assert 'TCGA' not in aliquot['gid'], 'Aliquot gid should be a uuid not {}'.format(aliquot['gid'])
            assert 'TCGA' in aliquot['data']['gdc_attributes']['submitter_id'], 'Aliquot gid should be a uuid not {}'.format(aliquot['gid'])


def test_simple(helpers, emitter_path_prefix):
    """ limit the result to a single project"""
    # see https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#example_2
    parameters = {
        "filters": {
            "op": "in",
            "content": {
                "field": "cases.submitter_id", "value": ["TCGA-02-0003"]
            }
        }

    }
    validate(helpers, emitter_path_prefix, parameters)
