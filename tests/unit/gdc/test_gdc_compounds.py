import os
import contextlib
import shutil
import pytest

from transform.gdc.compounds import transform


@pytest.fixture
def compounds_path(request):
    """get the full path of the test input"""
    return os.path.join(request.fspath.dirname, 'source/gdc/compounds/*.tsv')


@pytest.fixture
def id_lookup_path(request):
    """get the full path of the test input"""
    return os.path.join(request.fspath.dirname, 'source/gdc/id_lookup.tsv')


@pytest.fixture
def project_lookup_path(request):
    """get the full path of the test input"""
    return os.path.join(request.fspath.dirname, 'source/gdc/project_lookup.tsv')


def validate(helpers, emitter_directory, compounds_path, id_lookup_path, project_lookup_path):
    """ run xform and test results"""
    compound_file = os.path.join(emitter_directory, 'Compound.Vertex.json.gz')

    compounds_edge_file = os.path.join(emitter_directory, 'Case_Compounds_Compound.Edge.json.gz')
    cases_edge_file = os.path.join(emitter_directory, 'Compound_Cases_Case.Edge.json.gz')
    cpp_edge_file = os.path.join(emitter_directory, 'Compound_Projects_Project.Edge.json.gz')
    pcc_edge_file = os.path.join(emitter_directory, 'Project_Compounds_Compound.Edge.json.gz')

    all_files = [
        # vertices
        compound_file,
        # edges
        cases_edge_file, compounds_edge_file, cpp_edge_file, pcc_edge_file
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(compounds=compounds_path,
              id_lookup_path=id_lookup_path,
              project_lookup_path=project_lookup_path,
              emitter_prefix=None,
              emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)


def test_simple(helpers, emitter_directory, compounds_path, id_lookup_path, project_lookup_path):

    validate(helpers, emitter_directory, compounds_path, id_lookup_path, project_lookup_path)
