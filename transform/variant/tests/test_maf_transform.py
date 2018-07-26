#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" test maf_transform """

import pytest
import maf_transform
import os
import contextlib
import json

import requests_cache

# cache responses
requests_cache.install_cache('tests/test.myvariant.info',
                             allowable_codes=(200, 400, 404),
                             allowable_methods=('GET', 'POST'))


@pytest.fixture
def maf_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test.maf')


@pytest.fixture
def gz_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'gz-test.maf.gz')


@pytest.fixture
def no_center_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test-no-center.maf')


@pytest.fixture
def LARP6_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test-LARP6.maf')


@pytest.fixture
def TMX3_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test-TMX3.maf')


@pytest.fixture
def TPTE_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test-TPTE.maf')


@pytest.fixture
def AXIN1_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test-AXIN1.maf')


@pytest.fixture
def RIMS1_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test-RIMS1.maf')


@pytest.fixture
def ZFP3_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test-ZFP3.maf')



@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


def validate(maf_file, emitter_path_prefix):
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    allelecall_file = '{}.AlleleCall.Edge.json'.format(emitter_path_prefix)
    callset_file = '{}.Callset.Vertex.json'.format(emitter_path_prefix)
    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(allele_file)
        os.remove(allelecall_file)
        os.remove(callset_file)
    # create output
    maf_transform.convert(maf_file, emitter_path_prefix)
    # test.Allele.Vertex.json
    error_message = 'maf_transform.convert({}, {}) should create {}' \
                    .format(maf_file, emitter_path_prefix, allele_file)
    assert os.path.isfile(allele_file),  error_message
    # test allele contents
    with open(allele_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allele = json.loads(line)
            # minimum protograph keys
            assert list(allele.keys()) == ['gid', 'label', 'data'], \
                'expected keys'
            # should not be emoty
            for k in allele.keys():
                assert allele[k], 'empty key %s' % k

            # mandatory keys
            required_keys = ['genome', 'chromosome', 'start', 'end',
                             'reference_bases', 'alternate_bases']
            for k in required_keys:
                assert allele['data'][k], 'empty key %s' % k

            # optional keys, if set should be non null
            optional_keys = ['annotations', 'myvariantinfo']
            for k in optional_keys:
                if k in allele['data']:
                    assert allele['data'][k], 'empty key %s' % k

    # test.Callset.Vertex.json
    error_message = 'maf_transform.convert({}, {}) should create {}' \
                    .format(maf_file, emitter_path_prefix, callset_file)
    assert os.path.isfile(callset_file),  error_message
    # test allele contents
    with open(callset_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            callset = json.loads(line)
            # mandatory keys
            required_keys = ['tumor_biosample_id', 'normal_biosample_id',
                             'call_method']
            for k in required_keys:
                assert callset['data'][k], 'empty key %s' % k

    # test.AlleleCall.Edge.json
    error_message = 'maf_transform.convert({}, {}) should create {}' \
                    .format(maf_file, emitter_path_prefix, allelecall_file)
    assert os.path.isfile(allelecall_file),  error_message
    # test allele contents
    with open(allelecall_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allelecall = json.loads(line)


def test_simple(maf_file, emitter_path_prefix):
    """ simple test """
    validate(maf_file, emitter_path_prefix)


def test_gz(gz_file, emitter_path_prefix):
    """ simple test """
    validate(gz_file, emitter_path_prefix)


def test_get_value():
    """ test default return"""
    assert maf_transform.get_value({'foo': 0}, 'bar', 1) == 1


def test_no_center(no_center_file, emitter_path_prefix):
    """ 'Center column' renamed """
    validate(no_center_file, emitter_path_prefix)


def test_larp6(LARP6_file, emitter_path_prefix):
    """ no alt """
    validate(LARP6_file, emitter_path_prefix)


def test_tpte(TPTE_file, emitter_path_prefix):
    """ no alt """
    validate(TPTE_file, emitter_path_prefix)


def test_AXIN1(AXIN1_file, emitter_path_prefix):
    """ no alt """
    validate(AXIN1_file, emitter_path_prefix)


def test_TMX3(TMX3_file, emitter_path_prefix):
    """ no alt """
    validate(TMX3_file, emitter_path_prefix)


def test_RIMS1(RIMS1_file, emitter_path_prefix):
    """ no alt """
    validate(RIMS1_file, emitter_path_prefix)


def test_ZFP3(ZFP3_file, emitter_path_prefix):
    """ no alt """
    validate(ZFP3_file, emitter_path_prefix)
