
import os
import pytest
from transform.ccle.drug_response import transform


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


@pytest.fixture
def sample_info_file(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'source/ccle/DepMap-2018q3-celllines.csv')


def test_simple(helpers, emitter_path_prefix, sample_info_file):
    """ limit the result to a single project"""
    transform()
