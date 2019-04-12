""" test maf_transform """

import pytest
from transform.ccle.ccle_maf_transform import transform
from transform.ccle.ccle_maf_transform import CCLE_EXTENSION_CALLSET_KEYS
from bmeg.vertex import Allele, Callset, Gene, Aliquot, Sample, Case, Project
from bmeg.maf.maf_transform import STANDARD_MAF_KEYS
from bmeg.ioutils import reader

import os
import contextlib
import json


@pytest.fixture
def maf_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/ccle/mafs/*/vep.maf')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


@pytest.fixture
def ccle_sample_path(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'outputs/ccle/Sample.Vertex.json.gz')


def validate(helpers, maf_file, emitter_path_prefix, ccle_sample_path):
    allele_file = os.path.join(emitter_path_prefix, 'maf.Allele.Vertex.json.gz')
    allelecall_file = os.path.join(emitter_path_prefix, 'maf.SomaticVariant.Edge.json.gz')
    callset_file = os.path.join(emitter_path_prefix, 'maf.Callset.Vertex.json.gz')
    allelein_file = os.path.join(emitter_path_prefix, 'maf.AlleleIn.Edge.json.gz')
    callsetfor_file = os.path.join(emitter_path_prefix, 'maf.HasCallset.Edge.json.gz')
    aliquot_file = os.path.join(emitter_path_prefix, 'maf.Aliquot.Vertex.json.gz')
    aliquotfor_file = os.path.join(emitter_path_prefix, 'maf.HasAliquot.Edge.json.gz')
    sample_file = os.path.join(emitter_path_prefix, 'maf.Sample.Vertex.json.gz')
    samplefor_file = os.path.join(emitter_path_prefix, 'maf.HasSample.Edge.json.gz')
    case_file = os.path.join(emitter_path_prefix, 'maf.Case.Vertex.json.gz')
    project_file = os.path.join(emitter_path_prefix, 'maf.Project.Vertex.json.gz')
    inproject_file = os.path.join(emitter_path_prefix, 'maf.HasCase.Edge.json.gz')
    all_files = [allele_file, allelecall_file, callset_file, allelein_file, callsetfor_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)
    # create output
    transform(
        mafpath=maf_file,
        ccle_sample_path=ccle_sample_path,
        emitter_directory=emitter_path_prefix)

    # test/maf.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    # test/maf.Callset.Vertex.json
    callset_count = helpers.assert_vertex_file_valid(Callset, callset_file)
    # test/maf.AlleleIn.Edge.json
    helpers.assert_edge_file_valid(Allele, Gene, allelein_file)
    # test/maf.AlleleCall.Edge.json
    helpers.assert_edge_file_valid(Callset, Allele, allelecall_file)
    # test/maf.CallsetFor.Edge.json
    helpers.assert_edge_file_valid(Aliquot, Callset, callsetfor_file)
    # test/maf.Aliquot.Vertex.json
    helpers.assert_vertex_file_valid(Aliquot, aliquot_file)
    # test/maf.AliquotFor.Edge.json
    helpers.assert_edge_file_valid(Sample, Aliquot, aliquotfor_file)
    # test/maf.Sample.Vertex.json
    helpers.assert_vertex_file_valid(Sample, sample_file)
    # test/maf.SampleFor.Edge.json
    helpers.assert_edge_file_valid(Case, Sample, samplefor_file)
    # test/maf.Case.Vertex.json
    helpers.assert_vertex_file_valid(Case, case_file)
    # test/maf.Project.Vertex.json
    helpers.assert_vertex_file_valid(Project, project_file)
    # test/maf.InProject.Edge.json
    helpers.assert_edge_file_valid(Project, Case, inproject_file)

    assert callset_count > 0, 'There should be at least one callset'
    with reader(callset_file) as f:
        for line in f:
            # should be json
            callset = json.loads(line)
            # source should be ccle
            assert callset['data']['source'] == 'ccle', 'source should be ccle'
            # from & to should be ids, not gids
            assert 'Aliquot' not in callset['data']['tumor_aliquot_id'], 'tumor_aliquot_id should not have Aliquot gid'

    # test AlleleCall contents
    with reader(allelecall_file) as f:
        for line in f:
            # should be json
            allelecall = json.loads(line)
            # optional keys, if set should be non null
            for k in CCLE_EXTENSION_CALLSET_KEYS:
                if k in allelecall['data']:
                    assert allelecall['data'][k], 'empty key %s' % k

    # test Allele contents
    with reader(allele_file) as f:
        dbSNP_RS_count = 0
        for line in f:
            # should be json
            allele = json.loads(line)
            # ref & allele should be different
            assert allele['data']['reference_bases'] != allele['data']['alternate_bases'], 'reference should not equal alternate'
            for k in STANDARD_MAF_KEYS:
                if k in allele['data']:
                    assert allele['data'][k], 'empty key %s' % k
            if 'dbSNP_RS' in allele['data']:
                dbSNP_RS_count += 1
        assert dbSNP_RS_count > 0, 'should have some dbSNP_RS set'

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Gene', 'Aliquot']
    )


def test_simple(helpers, maf_file, emitter_path_prefix, ccle_sample_path):
    """ simple test """
    validate(helpers, maf_file, emitter_path_prefix, ccle_sample_path)
