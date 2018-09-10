
import os
import contextlib
import pytest
import json
from transform.gdc.cases import transform
from bmeg.vertex import Biosample, Aliquot, Individual, Project
from bmeg.emitter import JSONEmitter


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/')


def validate(helpers, emitter_path_prefix, parameters):
    """ run xform and test results"""
    biosample_file = '{}Biosample.Vertex.json'.format(emitter_path_prefix)
    aliquot_file = '{}Aliquot.Vertex.json'.format(emitter_path_prefix)
    individual_file = '{}Individual.Vertex.json'.format(emitter_path_prefix)

    biosamplefor_file = '{}BiosampleFor.Edge.json'.format(emitter_path_prefix)
    aliquotfor_file = '{}AliquotFor.Edge.json'.format(emitter_path_prefix)
    inproject_file = '{}InProject.Edge.json'.format(emitter_path_prefix)

    all_files = [biosample_file, aliquot_file, individual_file, biosamplefor_file, aliquotfor_file, inproject_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    emitter = JSONEmitter(emitter_path_prefix)
    transform(emitter, parameters)
    emitter.close()
    # test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Biosample, biosample_file)
    # test.Aliquot.Vertex.json
    helpers.assert_vertex_file_valid(Aliquot, aliquot_file)
    # test.Individual.Vertex.json
    helpers.assert_vertex_file_valid(Individual, individual_file)
    # test.BiosampleFor.Edge.json
    helpers.assert_edge_file_valid(Biosample, Individual, biosamplefor_file)
    # test.AliquotFor.Edge.json
    helpers.assert_edge_file_valid(Aliquot, Biosample, aliquotfor_file)
    # test.InProject.Edge.json
    helpers.assert_edge_file_valid(Individual, Project, inproject_file)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Project']
    )

    # test Aliquot contents
    with open(aliquot_file, 'r', encoding='utf-8') as f:
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
                "field": "submitter_id", "value": ["TCGA-02-0003"]
            }
        }

    }
    validate(helpers, emitter_path_prefix, parameters)
