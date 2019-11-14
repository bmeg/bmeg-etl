import contextlib
import os
import pytest
import shutil

from transform.pharmacodb.drug_response import transform


@pytest.fixture
def cells_path(request):
    return os.path.join(request.fspath.dirname, 'source/pharmacodb/cells.csv')


@pytest.fixture
def drugs_path(request):
    return os.path.join(request.fspath.dirname, 'source/pharmacodb/drugs.csv')


@pytest.fixture
def drug_annots_path(request):
    return os.path.join(request.fspath.dirname, 'source/pharmacodb/drug_annots.csv')


@pytest.fixture
def experiments_path(request):
    return os.path.join(request.fspath.dirname, 'source/pharmacodb/experiments.csv')


@pytest.fixture
def dose_response_path(request):
    return os.path.join(request.fspath.dirname, 'source/pharmacodb/dose_responses.csv')


@pytest.fixture
def profiles_path(request):
    return os.path.join(request.fspath.dirname, 'source/pharmacodb/profiles.csv')


@pytest.fixture
def cellline_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_id_lookup.tsv')


def validate(helpers, emitter_directory, cellline_lookup_path, cells_path, drugs_path,
             drug_annots_path, experiments_path, dose_response_path, profiles_path):

    drug_response_file = os.path.join(emitter_directory, 'DrugResponse.Vertex.json.gz')
    compound_file = os.path.join(emitter_directory, 'Compound.Vertex.json.gz')

    add_edge_file = os.path.join(emitter_directory, 'Aliquot_DrugResponse_DrugResponse.Edge.json.gz')
    daa_edge_file = os.path.join(emitter_directory, 'DrugResponse_Aliquot_Aliquot.Edge.json.gz')
    dcc_edge_file = os.path.join(emitter_directory, 'DrugResponse_Compounds_Compound.Edge.json.gz')
    cdd_edge_file = os.path.join(emitter_directory, 'Compound_DrugResponses_DrugResponse.Edge.json.gz')
    pcc_edge_file = os.path.join(emitter_directory, 'Project_Compounds_Compound.Edge.json.gz')
    cpp_edge_file = os.path.join(emitter_directory, 'Compound_Projects_Project.Edge.json.gz')

    all_files = [
        drug_response_file, compound_file,
        add_edge_file, daa_edge_file, dcc_edge_file, cdd_edge_file,
        pcc_edge_file, cpp_edge_file
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(
        cellline_lookup_path=cellline_lookup_path,
        cells_path=cells_path,
        drugs_path=drugs_path,
        drug_annots_path=drug_annots_path,
        experiments_path=experiments_path,
        dose_response_path=dose_response_path,
        profiles_path=profiles_path,
        emitter_prefix=None,
        emitter_directory=emitter_directory
    )

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot', 'Project'])


def test_simple(helpers, emitter_directory, cellline_lookup_path, cells_path, drugs_path,
                drug_annots_path, experiments_path, dose_response_path, profiles_path):

    validate(helpers, emitter_directory, cellline_lookup_path, cells_path, drugs_path,
             drug_annots_path, experiments_path, dose_response_path, profiles_path)
