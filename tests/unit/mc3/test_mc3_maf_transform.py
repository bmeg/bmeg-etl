

""" test maf_transform """

import pytest

import transform.mc3.mc3_maf_transform as mc3_maf_transform
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
    return os.path.join(request.fspath.dirname, 'source/mc3/tcga_test-NO_CENTER.maf')


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
    callset_file = os.path.join(emitter_directory, 'SomaticCallset.Vertex.json.gz')
    # deadletter_file = os.path.join(emitter_directory, 'Deadletter.Vertex.json.gz')

    aliquot_callset_edge_file = os.path.join(emitter_directory, 'Aliquot_SomaticCallsets_SomaticCallset.Edge.json.gz')
    callset_aliquot_edge_file = os.path.join(emitter_directory, 'SomaticCallset_Aliquots_Aliquot.Edge.json.gz')
    allele_callset_edge_file = os.path.join(emitter_directory, 'Allele_SomaticCallsets_SomaticCallset.Edge.json.gz')
    callset_allele_edge_file = os.path.join(emitter_directory, 'SomaticCallset_Alleles_Allele.Edge.json.gz')

    all_files = [allele_file, callset_file,
                 aliquot_callset_edge_file, callset_aliquot_edge_file,
                 allele_callset_edge_file, callset_allele_edge_file]

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

    # test alleles edge contents
    with reader(allele_callset_edge_file) as f:
        for line in f:
            # should be json
            allelecall = json.loads(line)
            if not (allelecall['from'].startswith("SomaticCallset") and allelecall['to'].startswith("Allele")):
                continue
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

    # check callset
    with reader(callset_file) as f:
        for line in f:
            # should be json
            callset = json.loads(line)
            assert callset['gid'].startswith('SomaticCallset:MC3:'), 'should start with SomaticCallset:MC3:xxx'
            assert not callset['gid'].startswith('SomaticCallset:MC3:Aliquot:'), 'should NOT start with SomaticCallset:MC3:Aliquot:xxx'
            assert callset['data']['tumor_aliquot_id'] != callset['data']['normal_aliquot_id'], 'tumor should not equal normal'
            assert 'Aliquot:' not in callset['data']['tumor_aliquot_id'], 'tumor_aliquot_id should not have Aliquot gid'
            assert 'Aliquot:' not in callset['data']['normal_aliquot_id'], 'normal_aliquot_id should not have Aliquot gid'

    # check callsetfor
    with reader(callset_aliquot_edge_file) as f:
        for line in f:
            # should be json
            callsetfor = json.loads(line)
            assert callsetfor['from'].startswith('SomaticCallset:MC3:'), 'from should be a callset'
            assert callsetfor['to'].startswith('Aliquot:'), 'to should be an aliquot'

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Aliquot']
    )
    return all_files


def test_simple(helpers, maf_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ simple test """
    validate(helpers, maf_file, emitter_directory, id_lookup_path, project_lookup_path)


def test_gz(helpers, gz_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ simple test """
    validate(helpers, gz_file, emitter_directory, id_lookup_path, project_lookup_path)


def test_no_center(helpers, no_center_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ 'Center column' renamed """
    validate(helpers, no_center_file, emitter_directory, id_lookup_path, project_lookup_path)


def test_NO_BARCODE(helpers, NO_BARCODE_file, emitter_directory, id_lookup_path, project_lookup_path):
    """ no barcode """
    validate(helpers, NO_BARCODE_file, emitter_directory, id_lookup_path, project_lookup_path)
