""" test maf_transform """

import os
import contextlib
import pytest
import logging
from transform.compound.transform import transform
from bmeg.vertex import Compound, DrugResponse

ALL_FILES = """
Compound.Vertex.json
ResponseTo.Edge.json
""".strip().split()


@pytest.fixture
def output_dir(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs')


def validate(helpers, output_dir, emitter_directory):

    all_files = [os.path.join(emitter_directory, f) for f in ALL_FILES]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)
    # create output
    transform(output_dir=output_dir, emitter_directory=emitter_directory)
    # check output
    compounds = all_files[0]
    response_tos = all_files[1]
    helpers.assert_vertex_file_valid(Compound, compounds)
    helpers.assert_edge_file_valid(DrugResponse, Compound, response_tos)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['DrugResponse'])


def test_simple(caplog, helpers, output_dir, emitter_directory):
    """ simple test """
    caplog.set_level(logging.INFO)
    validate(helpers, output_dir, emitter_directory)
