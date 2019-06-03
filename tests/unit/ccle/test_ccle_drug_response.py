
import os
import pytest
import contextlib
from transform.ccle.ccle_drug_response import transform
from bmeg.vertex import DrugResponse, Aliquot, Compound, Project


@pytest.fixture
def drug_response_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv')


@pytest.fixture
def cellline_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_lookup.tsv')


@pytest.fixture
def project_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_project_lookup.tsv')


def validate(helpers, emitter_directory, cellline_lookup_path, project_lookup_path, drug_response_path):
    """ run xform and test results"""
    profile_file = os.path.join(
        emitter_directory, 'drug_response.DrugResponse.Vertex.json.gz')
    profile_in_file = os.path.join(
        emitter_directory, 'drug_response.ResponseIn.Edge.json.gz')
    response_to_file = os.path.join(
        emitter_directory, 'drug_response.ResponseTo.Edge.json.gz')
    compound_file = os.path.join(
        emitter_directory, 'drug_response.Compound.Vertex.json.gz')
    tested_in_file = os.path.join(
        emitter_directory, 'drug_response.TestedIn.Edge.json.gz')

    all_files = [profile_file, profile_in_file, response_to_file, compound_file,
                 tested_in_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    transform(cellline_lookup_path=cellline_lookup_path,
              project_lookup_path=project_lookup_path,
              drug_response_path=drug_response_path,
              emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(DrugResponse, profile_file)
    helpers.assert_edge_file_valid(DrugResponse, Aliquot, profile_in_file)
    helpers.assert_edge_file_valid(DrugResponse, Compound, response_to_file)
    helpers.assert_edge_file_valid(Compound, Project, tested_in_file)
    helpers.assert_vertex_file_valid(Compound, compound_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot', 'Project'])


def test_simple(helpers, emitter_directory, cellline_lookup_path, project_lookup_path, drug_response_path):
    """ limit the result to a single project"""
    validate(helpers, emitter_directory, cellline_lookup_path, project_lookup_path, drug_response_path)
