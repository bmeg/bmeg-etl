

""" test maf_transform """

import pytest
import transform.mc3.mc3_maf_transform as mc3_maf_transform
from transform.mc3.mc3_maf_transform import MC3_EXTENSION_MAF_KEYS
from transform.mc3.mc3_maf_transform import MC3_EXTENSION_CALLSET_KEYS
from bmeg.maf.maf_transform import STANDARD_MAF_KEYS
from bmeg.maf.maf_transform import get_value
from bmeg.vertex import Allele, Callset, Gene, Aliquot
import os
import contextlib
import json


@pytest.fixture
def maf_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'tcga_test.maf')


@pytest.fixture
def gz_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'tcga_gz-test.maf.gz')


@pytest.fixture
def no_center_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'tcga_test-no-center.maf')


@pytest.fixture
def NO_REF_ALT_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'tcga_test-NO_REF_ALT.maf')


@pytest.fixture
def NO_BARCODE_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'tcga_test-NO_BARCODE.maf')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/test')


@pytest.fixture
def gdc_aliquot_path(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'outputs/gdc/gdc.Aliquot.Vertex.json')


def validate(helpers, maf_file, emitter_path_prefix, gdc_aliquot_path):
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    allelecall_file = '{}.AlleleCall.Edge.json'.format(emitter_path_prefix)
    callset_file = '{}.Callset.Vertex.json'.format(emitter_path_prefix)
    allelein_file = '{}.AlleleIn.Edge.json'.format(emitter_path_prefix)
    callsetfor_file = '{}.CallsetFor.Edge.json'.format(emitter_path_prefix)
    deadletter_file = '{}.Deadletter.Vertex.json'.format(emitter_path_prefix)
    all_files = [allele_file, allelecall_file, callset_file, allelein_file, callsetfor_file, deadletter_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)
    # create output
    mc3_maf_transform.transform(maf_file, emitter_path_prefix, gdc_aliquot_path=gdc_aliquot_path)

    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    # test.Callset.Vertex.json
    helpers.assert_vertex_file_valid(Callset, callset_file)
    # test/test.AlleleCall.Edge.json
    helpers.assert_edge_file_valid(Allele, Callset, allelecall_file)
    # test/test.AlleleIn.Edge.json
    helpers.assert_edge_file_valid(Allele, Gene, allelein_file)
    # test/test.CallsetFor.Edge.json
    helpers.assert_edge_file_valid(Callset, Aliquot, callsetfor_file)

    with open(callset_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            callset = json.loads(line)
            # source should be ccle
            assert callset['data']['source'] == 'mc3', 'source should be ccle'
            assert 'Aliquot:' not in callset['data']['tumor_aliquot_id'], 'tumor_aliquot_id should not have Aliquot gid'
            assert 'Aliquot:' not in callset['data']['normal_aliquot_id'], 'normal_aliquot_id should not have Aliquot gid'


    # test AlleleCall contents
    with open(allelecall_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allelecall = json.loads(line)
            # optional keys, if set should be non null
            for k in MC3_EXTENSION_CALLSET_KEYS:
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
            for k in MC3_EXTENSION_MAF_KEYS:
                if k in allele['data']['annotations']['mc3']:
                    assert allele['data']['annotations']['mc3'][k], 'empty key %s' % k

    # check callset
    with open(callset_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            callset = json.loads(line)
            assert callset['gid'].startswith('Callset:mc3:'), 'should start with Callset:mc3:xxx'
            assert not callset['gid'].startswith('Callset:mc3:Aliquot:'), 'should start with Callset:mc3:xxx'
            assert callset['data']['tumor_aliquot_id'] != callset['data']['normal_aliquot_id'], 'tumor should not equal normal'

    # check callsetfor
    with open(callsetfor_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            callsetfor = json.loads(line)
            assert callsetfor['from'].startswith('Callset:mc3:'), 'from should be a callset'
            assert callsetfor['to'].startswith('Aliquot:'), 'to should be an aliquot'

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Gene', 'Aliquot']
    )


def test_simple(helpers, maf_file, emitter_path_prefix, gdc_aliquot_path):
    """ simple test """
    validate(helpers, maf_file, emitter_path_prefix, gdc_aliquot_path)


def test_gz(helpers, gz_file, emitter_path_prefix, gdc_aliquot_path):
    """ simple test """
    validate(helpers, gz_file, emitter_path_prefix, gdc_aliquot_path)


def test_get_value():
    """ test default return"""
    assert get_value({'foo': 0}, 'bar', 1) == 1


def test_no_center(helpers, no_center_file, emitter_path_prefix, gdc_aliquot_path):
    """ 'Center column' renamed """
    validate(helpers, no_center_file, emitter_path_prefix, gdc_aliquot_path)


def test_NO_REF_ALT(helpers, NO_REF_ALT_file, emitter_path_prefix, gdc_aliquot_path):
    """ no start """
    with pytest.raises(AssertionError):
        validate(helpers, NO_REF_ALT_file, emitter_path_prefix, gdc_aliquot_path)
    deadletter_file = '{}.Deadletter.Vertex.json'.format(emitter_path_prefix)
    with open(deadletter_file, 'r', encoding='utf-8') as f:
        c = 0
        for line in f:
            json.loads(line)
            c += 1
        assert c == 1, 'We should have 1 dead letter'


def test_NO_BARCODE(helpers, NO_BARCODE_file, emitter_path_prefix, gdc_aliquot_path):
    """ no start """
    with pytest.raises(AssertionError):
        validate(helpers, NO_BARCODE_file, emitter_path_prefix, gdc_aliquot_path)
    deadletter_file = '{}.Deadletter.Vertex.json'.format(emitter_path_prefix)
    with open(deadletter_file, 'r', encoding='utf-8') as f:
        c = 0
        for line in f:
            json.loads(line)
            c += 1
        assert c == 1, 'We should have 1 dead letter'
