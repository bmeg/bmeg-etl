""" test maf_transform """

import pytest
from transform.ccle.ccle_gen3_file_transform import transform
from bmeg.vertex import File, Aliquot

import os
import contextlib


@pytest.fixture
def maf_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/ccle/mafs/*/vep.maf')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


@pytest.fixture
def ccle_biosample_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs/ccle/Biosample.Vertex.json.gz')


@pytest.fixture
def dvc_file(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs.ccle.maf.dvc')


def validate(helpers, maf_file, emitter_path_prefix, ccle_biosample_path, dvc_file):
    file = os.path.join(emitter_path_prefix, 'File.Vertex.json.gz')
    derived_from = os.path.join(emitter_path_prefix, 'DerivedFrom.Edge.json.gz')
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in [file, derived_from]:
            os.remove(f)
    # create output
    transform(
        mafpath=maf_file,
        ccle_biosample_path=ccle_biosample_path,
        emitter_directory=emitter_path_prefix,
        dvc_file=dvc_file)
    #
    helpers.assert_vertex_file_valid(File, file)
    helpers.assert_edge_file_valid(File, Aliquot, derived_from)


def test_simple(helpers, maf_file, emitter_path_prefix, ccle_biosample_path, dvc_file):
    """ simple test """
    validate(helpers, maf_file, emitter_path_prefix, ccle_biosample_path, dvc_file)
