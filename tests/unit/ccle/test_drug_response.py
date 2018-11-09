
import os
import pytest
import contextlib
from transform.ccle.drug_response import transform
from bmeg.vertex import ParamacalogicalProfile, Aliquot, Compound, Biosample, Individual, Project


@pytest.fixture
def drug_response_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv')


@pytest.fixture
def biosample_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs/ccle/Biosample.Vertex.json.gz')


def validate(helpers, emitter_directory, biosample_path, drug_response_path):
    """ run xform and test results"""
    profile_file = os.path.join(emitter_directory, 'drug_response.ParamacalogicalProfile.Vertex.json.gz')
    profile_in_file = os.path.join(emitter_directory, 'drug_response.ParamacalogicalProfileIn.Edge.json.gz')
    response_to_file = os.path.join(emitter_directory, 'drug_response.ResponseTo.Edge.json.gz')
    compound_file = os.path.join(emitter_directory, 'drug_response.Compound.Vertex.json.gz')

    all_files = [profile_file, profile_in_file, response_to_file, compound_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    transform(biosample_path=biosample_path,
              drug_response_path=drug_response_path,
              emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(ParamacalogicalProfile, profile_file)
    helpers.assert_edge_file_valid(ParamacalogicalProfile, Aliquot, profile_in_file)
    helpers.assert_edge_file_valid(ParamacalogicalProfile, Compound, response_to_file)
    helpers.assert_vertex_file_valid(Compound, compound_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])
    # missing vertexes
    for f in ['drug_response.Aliquot.Vertex.json.gz', 'drug_response.Biosample.Vertex.json.gz',
              'drug_response.Individual.Vertex.json.gz', 'drug_response.Project.Vertex.json.gz']:
        v = eval(f.split('.')[1])
        f = os.path.join(emitter_directory, f)
        helpers.assert_vertex_file_valid(v, f)
    # missing edges
    for f, v1, v2 in [('drug_response.AliquotFor.Edge.json.gz', Aliquot, Biosample),
                      ('drug_response.BiosampleFor.Edge.json.gz', Biosample, Individual),
                      ('drug_response.InProject.Edge.json.gz', Individual, Project)]:
        f = os.path.join(emitter_directory, f)
        helpers.assert_edge_file_valid(v1, v2, f)


def test_simple(helpers, emitter_directory, biosample_path, drug_response_path):
    """ limit the result to a single project"""
    validate(helpers, emitter_directory, biosample_path, drug_response_path)
