""" test maf_transform """

import os
import contextlib
import pytest
import logging
import shutil
from transform.compound.transform import transform
from bmeg.stores import new_store

ALL_FILES = """
normalized.Compound.Vertex.json.gz
normalized.DrugResponse_Compounds_Compound.Edge.json.gz
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
        shutil.rmtree(emitter_directory)
        os.remove(store_path)

    # create output
    transform(output_dir=output_dir,
              emitter_directory=emitter_directory)

    # check output
    compounds = all_files[0]
    response_tos = all_files[1]
    helpers.assert_vertex_file_valid(compounds)
    helpers.assert_edge_file_valid(response_tos)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['DrugResponse'])


def test_simple(caplog, helpers, output_dir, emitter_directory, store_path):
    """ simple test """
    caplog.set_level(logging.INFO)
    validate(helpers, output_dir, emitter_directory, store_path)
