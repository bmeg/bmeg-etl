""" test maf_transform """

import pytest
import transform.ccle.ccle_maf_transform as ccle_maf_transform
from transform.ccle.ccle_maf_transform import CCLE_EXTENSION_CALLSET_KEYS, CCLE_EXTENSION_MAF_KEYS
from bmeg.vertex import Allele, Callset, Gene, Aliquot
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
    return os.path.join(request.fspath.dirname, 'test/test')


def validate(helpers, maf_file, emitter_path_prefix, harvest=True, filter=[]):
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    allelecall_file = '{}.AlleleCall.Edge.json'.format(emitter_path_prefix)
    callset_file = '{}.Callset.Vertex.json'.format(emitter_path_prefix)
    allelein_file = '{}.AlleleIn.Edge.json'.format(emitter_path_prefix)
    callsetfor_file = '{}.CallsetFor.Edge.json'.format(emitter_path_prefix)
    all_files = [allele_file, allelecall_file, callset_file, allelein_file, callsetfor_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)
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
    # test/test.CallsetFor.Edge.json
    helpers.assert_edge_file_valid(Callset, Aliquot, callsetfor_file)

    with open(callset_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            callset = json.loads(line)
            # source should be ccle
            assert callset['data']['source'] == 'ccle', 'source should be ccle'
            # from & to should be ids, not gids
            assert 'Aliquot' not in callset['data']['tumor_aliquot_id'], 'tumor_aliquot_id should not have Aliquot gid'

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
            # ref & allele should be different
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
        all_files,
        exclude_labels=['Gene', 'Aliquot']
    )


def test_simple(helpers, maf_file, emitter_path_prefix):
    """ simple test """
    validate(helpers, maf_file, emitter_path_prefix)
