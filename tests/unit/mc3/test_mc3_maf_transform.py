

""" test maf_transform """

import pytest

import transform.mc3.mc3_maf_transform as mc3_maf_transform
from transform.mc3.mc3_maf_transform import MC3_EXTENSION_CALLSET_KEYS
from bmeg.maf.maf_transform import STANDARD_MAF_KEYS
from bmeg.maf.maf_transform import get_value
from bmeg.ioutils import reader

import os
import contextlib
import shutil
import json


@pytest.fixture
def maf_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/mc3/tcga_test.maf')


@pytest.fixture
def gz_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/mc3/tcga_gz-test.maf.gz')


@pytest.fixture
def no_center_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/mc3/tcga_test-no-center.maf')


@pytest.fixture
def NO_REF_ALT_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/mc3/tcga_test-NO_REF_ALT.maf')


@pytest.fixture
def NO_BARCODE_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/mc3/tcga_test-NO_BARCODE.maf')


@pytest.fixture
def id_lookup_path(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdc/id_lookup.tsv')


@pytest.fixture
def project_lookup_path(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdc/project_lookup.tsv')


def validate(helpers, maf_file, emitter_directory, id_lookup_path, project_lookup_path):
    allele_file = os.path.join(emitter_directory, 'Allele.Vertex.json.gz')
    callset_file = os.path.join(emitter_directory, 'Callset.Vertex.json.gz')
    # deadletter_file = os.path.join(emitter_directory, 'Deadletter.Vertex.json.gz')

    callsets_edge_file = os.path.join(emitter_directory, 'callsets.Edge.json.gz')
    aliquots_edge_file = os.path.join(emitter_directory, 'aliquots.Edge.json.gz')
    alleles_edge_file = os.path.join(emitter_directory, 'alleles.Edge.json.gz')

    all_files = [
        allele_file, callset_file,  # deadletter_file,
        callsets_edge_file, aliquots_edge_file, alleles_edge_file
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    mc3_maf_transform.transform(mafpath=maf_file,
                                id_lookup_path=id_lookup_path,
                                project_lookup_path=project_lookup_path,
                                emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Gene', 'Aliquot']
    )

    # test alleles edge contents
    with reader(alleles_edge_file) as f:
        for line in f:
            # should be json
            allelecall = json.loads(line)
            if not (allelecall['from'].startswith("Callset") and allelecall['to'].startswith("Allele")):
                continue
            # optional keys, if set should be non null
            for k in MC3_EXTENSION_CALLSET_KEYS:
                if k in allelecall['data']:
                    assert allelecall['data'][k], 'empty key %s' % k
            assert '|' not in allelecall['data']['methods'], 'call_method should not have a | separator'
            allelecall_methods = set(allelecall['data']['methods'])
            possible_allelecall_methods = set(["RADIA", "MUTECT", "MUSE", "VARSCANS", "INDELOCATOR", "VARSCANI", "PINDEL", "SOMATICSNIPER"])
            assert allelecall_methods < possible_allelecall_methods, 'call_method should belong to vocabulary'

    # test Allele contents
    with reader(allele_file) as f:
        for line in f:
            # should be json
            allele = json.loads(line)
            assert allele['data']['reference_bases'] != allele['data']['alternate_bases'], 'reference should not equal alternate'
            for k in STANDARD_MAF_KEYS:
                if k in allele['data']:
                    assert allele['data'][k], 'empty key %s' % k

    # check callset
    with reader(callset_file) as f:
        for line in f:
            # should be json
            callset = json.loads(line)
            assert callset['gid'].startswith('Callset:MC3:'), 'should start with Callset:MC3:xxx'
            assert not callset['gid'].startswith('Callset:MC3:Aliquot:'), 'should NOT start with Callset:MC3:Aliquot:xxx'
            assert callset['data']['tumor_aliquot_id'] != callset['data']['normal_aliquot_id'], 'tumor should not equal normal'
            assert 'Aliquot:' not in callset['data']['tumor_aliquot_id'], 'tumor_aliquot_id should not have Aliquot gid'
            assert 'Aliquot:' not in callset['data']['normal_aliquot_id'], 'normal_aliquot_id should not have Aliquot gid'

    # check callsetfor
    with reader(aliquots_edge_file) as f:
        for line in f:
            # should be json
            callsetfor = json.loads(line)
            assert callsetfor['from'].startswith('Callset:MC3:'), 'from should be a callset'
            assert callsetfor['to'].startswith('Aliquot:'), 'to should be an aliquot'

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Gene', 'Aliquot']
    )
    return all_files


def test_simple(helpers, maf_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ simple test """
    validate(helpers, maf_file, emitter_directory, id_lookup_path, project_lookup_path)


def test_gz(helpers, gz_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ simple test """
    validate(helpers, gz_file, emitter_directory, id_lookup_path, project_lookup_path)


def test_get_value():
    """ test default return"""
    assert get_value({'foo': 0}, 'bar', 1) == 1


def test_no_center(helpers, no_center_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ 'Center column' renamed """
    validate(helpers, no_center_file, emitter_directory, id_lookup_path, project_lookup_path)


def test_NO_REF_ALT(helpers, NO_REF_ALT_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ no start """
    with pytest.raises(AssertionError):
        validate(helpers, NO_REF_ALT_file, emitter_directory, id_lookup_path, project_lookup_path)
    deadletter_file = os.path.join(emitter_directory, 'Deadletter.Vertex.json.gz')
    with reader(deadletter_file) as f:
        c = 0
        for line in f:
            json.loads(line)
            c += 1
        assert c == 1, 'We should have 1 dead letter'


def test_NO_BARCODE(helpers, NO_BARCODE_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ no start """
    with pytest.raises(AssertionError):
        validate(helpers, NO_BARCODE_file, emitter_directory, id_lookup_path, project_lookup_path)
    deadletter_file = os.path.join(emitter_directory, 'Deadletter.Vertex.json.gz')
    with reader(deadletter_file) as f:
        c = 0
        for line in f:
            json.loads(line)
            c += 1
        assert c == 1, 'We should have 1 dead letter'
