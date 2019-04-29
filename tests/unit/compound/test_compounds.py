""" test maf_transform """

import os
import contextlib
import pytest
import logging
from transform.compound.transform import transform
from bmeg.vertex import Compound, DrugResponse
from bmeg.stores import new_store

ALL_FILES = """
normalized.Compound.Vertex.json.gz
normalized.ResponseTo.Edge.json.gz
""".strip().split()

VERTEX_FILES = """
test.Compound.Vertex.json
""".strip().split()

EDGE_FILES = """
test.DrugResponseIn.Edge.json
test.ResponseTo.Edge.json
""".strip().split()


@pytest.fixture
def output_dir(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs')


@pytest.fixture
def store_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/compound/sqlite.db')


def validate(helpers, output_dir, emitter_directory, store_path):

    all_files = [os.path.join(emitter_directory, f) for f in ALL_FILES]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)
        os.remove(store_path)
    # create output
    vertex_files = [os.path.join(output_dir, f) for f in VERTEX_FILES]
    edge_files = [os.path.join(output_dir, f) for f in EDGE_FILES]
    transform(vertex_files=vertex_files, edge_files=edge_files, output_dir=output_dir, emitter_directory=emitter_directory, store_path=store_path)
    # check output
    compounds = all_files[0]
    response_tos = all_files[1]
    helpers.assert_vertex_file_valid(Compound, compounds)
    helpers.assert_edge_file_valid(DrugResponse, Compound, response_tos)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['DrugResponse'])
    # ensure the store was created
    store = new_store('key-val', path=store_path)
    assert len([c for c in store.all()]) == 9, 'store should have 9 names'


def test_simple(caplog, helpers, output_dir, emitter_directory, store_path):
    """ simple test """
    caplog.set_level(logging.INFO)
    validate(helpers, output_dir, emitter_directory, store_path)
