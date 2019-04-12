
import os
import pytest
import contextlib
from transform.ccle.drug_response import transform
from bmeg.vertex import DrugResponse, Aliquot, Compound, Sample, Case, Project, Program


@pytest.fixture
def drug_response_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv')


@pytest.fixture
def sample_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs/ccle/Sample.Vertex.json.gz')


def validate(helpers, emitter_directory, sample_path, drug_response_path):
    """ run xform and test results"""
    profile_file = os.path.join(
        emitter_directory, 'drug_response.DrugResponse.Vertex.json.gz')
    profile_in_file = os.path.join(
        emitter_directory, 'drug_response.ResponseIn.Edge.json.gz')
    response_to_file = os.path.join(
        emitter_directory, 'drug_response.ResponseTo.Edge.json.gz')
    compound_file = os.path.join(
        emitter_directory, 'drug_response.Compound.Vertex.json.gz')

    all_files = [profile_file, profile_in_file,
                 response_to_file, compound_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    transform(sample_path=sample_path,
              drug_response_path=drug_response_path,
              emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(DrugResponse, profile_file)
    helpers.assert_edge_file_valid(DrugResponse, Aliquot, profile_in_file)
    helpers.assert_edge_file_valid(DrugResponse, Compound, response_to_file)
    helpers.assert_vertex_file_valid(Compound, compound_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])
    # missing vertexes
    for f in ['drug_response.Aliquot.Vertex.json.gz', 'drug_response.Sample.Vertex.json.gz',
              'drug_response.Case.Vertex.json.gz', 'drug_response.Project.Vertex.json.gz',
              'drug_response.Program.Vertex.json.gz']:
        v = eval(f.split('.')[1])
        f = os.path.join(emitter_directory, f)
        helpers.assert_vertex_file_valid(v, f)
    # missing edges
    for f, v1, v2 in [('drug_response.HasAliquot.Edge.json.gz', Sample, Aliquot),
                      ('drug_response.HasSample.Edge.json.gz', Case, Sample),
                      ('drug_response.HasCase.Edge.json.gz', Project, Case),
                      ('drug_response.HasProject.Edge.json.gz', Program, Project)]:
        f = os.path.join(emitter_directory, f)
        helpers.assert_edge_file_valid(v1, v2, f)


def test_simple(helpers, emitter_directory, sample_path, drug_response_path):
    """ limit the result to a single project"""
    validate(helpers, emitter_directory, sample_path, drug_response_path)
