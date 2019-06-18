import os
import shutil
import pytest

from transform.gdc.compounds import transform


@pytest.fixture
def compounds_path(request):
    """get the full path of the test input"""
    return os.path.join(request.fspath.dirname, 'source/gdc/compounds/*.tsv')


def validate(helpers, emitter_directory, compounds_path):
    """ run xform and test results"""
    compound_file = os.path.join(emitter_directory, 'Compound.Vertex.json.gz')

    compounds_edge_file = os.path.join(emitter_directory, 'compounds.Edge.json.gz')
    cases_edge_file = os.path.join(emitter_directory, 'cases.Edge.json.gz')
    projects_edge_file = os.path.join(emitter_directory, 'projects.Edge.json.gz')

    all_files = [
        # vertices
        compound_file,
        # edges
        projects_edge_file, cases_edge_file, compounds_edge_file
    ]

    # remove output
    shutil.rmtree(emitter_directory)

    # create output
    transform(input_path=compounds_path,
              emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)


def test_simple(helpers, emitter_directory, compounds_path):

    validate(helpers, emitter_directory, compounds_path)
