""" test maf_transform """

import os
import contextlib
import pytest
import logging
from transform.allele.transform import transform
from bmeg.vertex import Allele


@pytest.fixture
def output_directory(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'outputs')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/test')


def validate(helpers, output_directory, emitter_path_prefix):
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(allele_file)

    # create output
    transform(output_directory, prefix=emitter_path_prefix)
    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)


def test_simple(caplog, helpers, output_directory, emitter_path_prefix):
    """ simple test """
    caplog.set_level(logging.DEBUG)
    validate(helpers, output_directory, emitter_path_prefix)
