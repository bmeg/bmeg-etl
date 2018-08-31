""" test maf_transform """

import os
import contextlib
import pytest
import logging
import json
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


@pytest.fixture
def myvariantinfo_path(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/myvariantinfo/biothings_current_old_hg19.json.gz')


def validate_myvariantinfo_count(allele_file):
    myvariantinfo_count = 0
    with open(allele_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allele = json.loads(line)
            # optional keys, if set should be non null
            if allele['data']['annotations']['myvariantinfo']:
                myvariantinfo_count += 1
    assert myvariantinfo_count == 45, 'we should have had myvariantinfo hits'


def validate(helpers, output_directory, emitter_path_prefix, myvariantinfo_path):
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(allele_file)

    # check using memory store
    transform(output_directory,
              prefix=emitter_path_prefix,
              myvariant_store_name='myvariantinfo-memory',
              myvariantinfo_path=myvariantinfo_path)
    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    validate_myvariantinfo_count(allele_file)

    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(allele_file)

    # check using sqllite
    transform(output_directory,
              prefix=emitter_path_prefix,
              myvariant_store_name='myvariantinfo-sqlite',
              myvariantinfo_path=myvariantinfo_path)
    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    validate_myvariantinfo_count(allele_file)


def test_simple(caplog, helpers, output_directory, emitter_path_prefix, myvariantinfo_path):
    """ simple test """
    caplog.set_level(logging.DEBUG)
    validate(helpers, output_directory, emitter_path_prefix, myvariantinfo_path)
