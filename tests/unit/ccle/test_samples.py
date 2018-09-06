
import os
import contextlib
import pytest
from transform.ccle.samples import transform
from bmeg.vertex import Biosample, Aliquot


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/test')


@pytest.fixture
def sample_info_file(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_sample_info_file_2012-10-18.txt')


def validate(helpers, emitter_path_prefix, sample_info_file):
    """ run xform and test results"""
    biosample_file = '{}.Biosample.Vertex.json'.format(emitter_path_prefix)
    aliquot_file = '{}.Aliquot.Vertex.json'.format(emitter_path_prefix)
    aliquotfor_file = '{}.AliquotFor.Edge.json'.format(emitter_path_prefix)

    all_files = [biosample_file, aliquot_file, aliquotfor_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform(path=sample_info_file, prefix=emitter_path_prefix)

    # test.Biosample.Vertex.json
    helpers.assert_vertex_file_valid(Biosample, biosample_file)
    # test.Aliquot.Vertex.json
    helpers.assert_vertex_file_valid(Aliquot, aliquot_file)
    # test.AliquotFor.Edge.json
    helpers.assert_edge_file_valid(Aliquot, Biosample, aliquotfor_file)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(all_files)


def test_simple(helpers, emitter_path_prefix, sample_info_file):
    """ limit the result to a single project"""
    validate(helpers, emitter_path_prefix, sample_info_file)
