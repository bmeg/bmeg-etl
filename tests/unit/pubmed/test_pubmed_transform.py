
""" test maf_transform """

import os
import contextlib
import pytest
import json

from bmeg.emitter import JSONEmitter
from transform.pubmed.pubmed import parse_pubmed
from bmeg.vertex import Publication


@pytest.fixture
def pubmed_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'pubmed_test.xml')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


def validate(pubmed_file, emitter_path_prefix):
    publication_file = os.path.join(emitter_path_prefix, 'Publication.Vertex.json')
    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(publication_file)
    emitter = JSONEmitter(emitter_path_prefix)
    with open(pubmed_file) as handle:
        parse_pubmed(handle, emitter)
    assert_vertex_file_valid(Publication, publication_file)


def test_simple(pubmed_file, emitter_path_prefix):
    validate(pubmed_file, emitter_path_prefix)


# TODO - integrate test helpers ( currently on g2p branch)
def assert_data_keys_populated(data_class, vertex_dict):
    """ ensure that all non Union(NoneType,...) fields are not empty. """
    # mandatory keys
    for k in data_class.__dataclass_fields__.keys():
        field = data_class.__dataclass_fields__[k]
        # skip if union(None, ...)
        if 'typing.Union' in str(field.type) and 'NoneType' in str(field.type):
            continue
        assert vertex_dict['data'][k], 'empty key %s' % k


def assert_vertex_keys_populated(vertex_dict):
    """ ensure that graph keys populated """
    # minimum graph keys
    assert list(vertex_dict.keys()) == ['_id', 'gid', 'label', 'data'], \
        'expected keys'
    for k in vertex_dict.keys():
        assert vertex_dict[k], 'empty key %s' % k


def assert_vertex_file_valid(data_class, vertex_file_path):
    """ ensure file exists; populated with json objs with required keys """
    error_message = 'data_class {} {} ' \
                    .format(data_class, vertex_file_path)
    assert os.path.isfile(vertex_file_path), error_message
    with open(vertex_file_path, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            vertex_dict = json.loads(line)
            # should have all vertex keys
            assert_vertex_keys_populated(vertex_dict)
            assert_data_keys_populated(data_class, vertex_dict)
