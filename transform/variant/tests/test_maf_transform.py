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
def NO_REF_ALT_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'test-NO_REF_ALT.maf')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


def validate(maf_file, emitter_path_prefix, harvest=True, filter=[]):
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
    maf_transform.convert(maf_file, emitter_path_prefix,
                          harvest=harvest, filter=filter)
    # test.Allele.Vertex.json
    error_message = 'maf_transform.convert({}, {}) should create {}' \
                    .format(maf_file, emitter_path_prefix, allele_file)
    assert os.path.isfile(allele_file),  error_message
    # test allele contents
    with open(allele_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allele = json.loads(line)
            # minimum graph keys
            assert list(allele.keys()) == ['gid', 'label', 'data'], \
                'expected keys'
            # should not be empty
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
            # optional keys, if set should be non null
            optional_keys = ['t_depth', 't_ref_count', 't_alt_count',
                             'n_depth', 'n_ref_count', 'n_alt_count', 'FILTER',
                             'Match_Norm_Seq_Allele1',
                             'Match_Norm_Seq_Allele2',
                             'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
            for k in optional_keys:
                if k in allelecall['data']['info']:
                    assert allelecall['data']['info'][k], 'empty key %s' % k
    # test allelein (gene) contents
    with open(allelein_file, 'r', encoding='utf-8') as f:
        for line in f:
            # should be json
            allelein = json.loads(line)
            # minimum graph keys
            assert list(allelein.keys()) == \
                ['gid', 'label', 'from', 'to', 'data'], 'expected keys'


def test_simple(maf_file, emitter_path_prefix):
    """ simple test """
    validate(maf_file, emitter_path_prefix)


def test_simple_harvest(maf_file, emitter_path_prefix):
    """ simple test """
    validate(maf_file, emitter_path_prefix, harvest=False)


def test_simple_filter(maf_file, emitter_path_prefix):
    """ simple test """
    validate(maf_file, emitter_path_prefix,
             filter=['Allele:733057c6dcfcdaecc00a0212b65ee0421737c36f'])


def test_gz(gz_file, emitter_path_prefix):
    """ simple test """
    validate(gz_file, emitter_path_prefix)


def test_get_value():
    """ test default return"""
    assert maf_transform.get_value({'foo': 0}, 'bar', 1) == 1


def test_no_center(no_center_file, emitter_path_prefix):
    """ 'Center column' renamed """
    validate(no_center_file, emitter_path_prefix)


def test_NO_REF_ALT(NO_REF_ALT_file, emitter_path_prefix):
    """ no start """
    with pytest.raises(ValueError):
        validate(NO_REF_ALT_file, emitter_path_prefix)
