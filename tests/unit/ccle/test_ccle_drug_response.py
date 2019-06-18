
import os
import pytest
import shutil
from transform.ccle.ccle_drug_response import transform


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
    drug_response_file = os.path.join(emitter_directory, 'drug_response.DrugResponse.Vertex.json.gz')
    compound_file = os.path.join(emitter_directory, 'drug_response.Compound.Vertex.json.gz')

    drug_responses_edge_file = os.path.join(emitter_directory, 'drug_response.drug_responses.Edge.json.gz')
    drug_response_edge_file = os.path.join(emitter_directory, 'drug_response.drug_response.Edge.json.gz')
    compounds_edge_file = os.path.join(emitter_directory, 'drug_response.compounds.Edge.json.gz')
    projects_edge_file = os.path.join(emitter_directory, 'drug_response.projects.Edge.json.gz')

    all_files = [
        drug_response_file, compound_file,
        drug_responses_edge_file, drug_response_edge_file,
        compounds_edge_file, projects_edge_file
    ]

    # remove output
    shutil.rmtree(emitter_directory)

    transform(cellline_lookup_path=cellline_lookup_path,
              project_lookup_path=project_lookup_path,
              drug_response_path=drug_response_path,
              emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot', 'Project'])


def test_simple(helpers, emitter_directory, cellline_lookup_path, project_lookup_path, drug_response_path):
    """ limit the result to a single project"""
    validate(helpers, emitter_directory, cellline_lookup_path, project_lookup_path, drug_response_path)
