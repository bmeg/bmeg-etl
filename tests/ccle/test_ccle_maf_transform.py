""" test maf_transform """

import pytest
import transform.ccle.ccle_maf_transform as ccle_maf_transform
from bmeg.vertex import Allele, Callset, Gene

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
            optional_keys = ['cDNA_Change',
                             'Codon_Change',
                             'Protein_Change',
                             'isDeleterious',
                             'isTCGAhotspot',
                             'TCGAhsCnt',
                             'isCOSMIChotspot',
                             'COSMIChsCnt',
                             'ExAC_AF',
                             'WES_AC',
                             'WGS_AC',
                             'SangerWES_AC',
                             'SangerRecalibWES_AC',
                             'RNAseq_AC',
                             'HC_AC',
                             'RD_AC', ]
            for k in optional_keys:
                if k in allelecall['data']['info']:
                    assert allelecall['data']['info'][k], 'empty key %s' % k
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        [allele_file, allelecall_file, callset_file, allelein_file],
        exclude_labels=['Gene']
    )


def test_simple(helpers, maf_file, emitter_path_prefix):
    """ simple test """
    validate(helpers, maf_file, emitter_path_prefix)
