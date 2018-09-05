""" test maf_transform """

import pytest
import transform.ccle.ccle_maf_transform as ccle_maf_transform
from transform.ccle.ccle_maf_transform import CCLE_EXTENSION_CALLSET_KEYS, CCLE_EXTENSION_MAF_KEYS
from bmeg.vertex import Allele, Callset, Gene
from bmeg.maf.maf_transform import STANDARD_MAF_KEYS

import os
import contextlib
import json


@pytest.fixture
def maf_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'ccle_test.maf')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


def validate(helpers, maf_file, emitter_path_prefix, harvest=True, filter=[]):
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    allelecall_file = '{}.AlleleCall.Edge.json'.format(emitter_path_prefix)
    callset_file = '{}.Callset.Vertex.json'.format(emitter_path_prefix)
    allelein_file = '{}.AlleleIn.Edge.json'.format(emitter_path_prefix)
    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(allele_file)
        os.remove(allelecall_file)
        os.remove(callset_file)
        os.remove(allelein_file)
    # create output
    ccle_maf_transform.transform(maf_file, emitter_path_prefix)

    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    # test.Callset.Vertex.json
    helpers.assert_vertex_file_valid(Callset, callset_file)
    # test/test.AlleleIn.Edge.json
    helpers.assert_edge_file_valid(Allele, Gene, allelein_file)
    # test/test.AlleleCall.Edge.json
    helpers.assert_edge_file_valid(Allele, Callset, allelecall_file)
    # test.AlleleCall.Edge.json
    error_message = 'maf_transform.convert({}, {}) should create {}' \
                    .format(maf_file, emitter_path_prefix, allelecall_file)
    assert os.path.isfile(allelecall_file), error_message
    # test AlleleCall contents
    with open(allelecall_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allelecall = json.loads(line)
            # optional keys, if set should be non null
            for k in CCLE_EXTENSION_CALLSET_KEYS:
                if k in allelecall['data']['info']:
                    assert allelecall['data']['info'][k], 'empty key %s' % k
    # test Allele contents
    with open(allele_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allele = json.loads(line)
            assert allele['data']['reference_bases'] != allele['data']['alternate_bases'], 'reference should not equal alternate'            
            for k in STANDARD_MAF_KEYS:
                if k in allele['data']['annotations']['maf']:
                    assert allele['data']['annotations']['maf'][k], 'empty key %s' % k
            # optional keys, if set should be non null
            for k in CCLE_EXTENSION_MAF_KEYS:
                if k in allele['data']['annotations']['ccle']:
                    assert allele['data']['annotations']['ccle'][k], 'empty key %s' % k

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        [allele_file, allelecall_file, callset_file, allelein_file],
        exclude_labels=['Gene']
    )


def test_simple(helpers, maf_file, emitter_path_prefix):
    """ simple test """
    validate(helpers, maf_file, emitter_path_prefix)
