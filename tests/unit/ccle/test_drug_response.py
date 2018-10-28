
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
    profile_file = os.path.join(emitter_directory, 'ParamacalogicalProfile.Vertex.json')
    profile_in_file = os.path.join(emitter_directory, 'ParamacalogicalProfileIn.Edge.json')
    response_to_file = os.path.join(emitter_directory, 'ResponseTo.Edge.json')
    all_files = [profile_file, profile_in_file, response_to_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    transform(biosample_path=biosample_path,
              drug_response_path=drug_response_path,
              emitter_prefix=None,
              emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(ParamacalogicalProfile, profile_file)
    helpers.assert_edge_file_valid(ParamacalogicalProfile, Aliquot, profile_in_file)
    helpers.assert_edge_file_valid(ParamacalogicalProfile, Compound, response_to_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot', 'Compound'])
    # missing vertexes
    for f in ['Aliquot.Vertex.json', 'Biosample.Vertex.json', 'Compound.Vertex.json', 'Individual.Vertex.json', 'Project.Vertex.json']:
        v = eval(f.split('.')[0])
        f = os.path.join(emitter_directory, f)
        helpers.assert_vertex_file_valid(v, f)
    # missing edges
    for f, v1, v2 in [('AliquotFor.Edge.json', Aliquot, Biosample),
                      ('BiosampleFor.Edge.json', Biosample, Individual),
                      ('InProject.Edge.json', Individual, Project)]:
        f = os.path.join(emitter_directory, f)
        helpers.assert_edge_file_valid(v1, v2, f)


def test_simple(helpers, emitter_directory, biosample_path, drug_response_path):
    """ limit the result to a single project"""
    validate(helpers, emitter_directory, biosample_path, drug_response_path)
