""" test maf_transform """

import pytest
from transform.ccle.ccle_maf_transform import transform
from bmeg.ioutils import reader

import os
import contextlib
import shutil
import json


@pytest.fixture
def maf_file(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/test.maf')


@pytest.fixture
def cellline_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_id_lookup.tsv')


def validate(helpers, emitter_directory, maf_file, cellline_lookup_path):
    allele_file = os.path.join(emitter_directory, 'maf.Allele.Vertex.json.gz')
    callset_file = os.path.join(emitter_directory, 'maf.SomaticCallset.Vertex.json.gz')

    aliquot_callset_edge_file = os.path.join(emitter_directory, 'maf.Aliquot_SomaticCallsets_SomaticCallset.Edge.json.gz')
    callset_aliquot_edge_file = os.path.join(emitter_directory, 'maf.SomaticCallset_Aliquots_Aliquot.Edge.json.gz')
    allele_callset_edge_file = os.path.join(emitter_directory, 'maf.Allele_SomaticCallsets_SomaticCallset.Edge.json.gz')
    callset_allele_edge_file = os.path.join(emitter_directory, 'maf.SomaticCallset_Alleles_Allele.Edge.json.gz')

    all_files = [allele_file, callset_file,
                 aliquot_callset_edge_file, callset_aliquot_edge_file,
                 allele_callset_edge_file, callset_allele_edge_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(
        mafpath=maf_file,
        cellline_lookup_path=cellline_lookup_path,
        emitter_directory=emitter_directory
    )

    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    with reader(callset_file) as f:
        for line in f:
            # should be json
            callset = json.loads(line)
            # source should be ccle
            assert "CCLE" in callset['gid'], 'gid should contain CCLE'
            # from & to should be ids, not gids
            assert 'Aliquot' not in callset['data']['tumor_aliquot_id'], 'tumor_aliquot_id should not have Aliquot gid'

    # test Allele contents
    with reader(allele_file) as f:
        for line in f:
            # should be json
            allele = json.loads(line)
            # ref & allele should be different
            assert allele['data']['reference_bases'] != allele['data']['alternate_bases'], 'reference should not equal alternate'

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Aliquot']
    )


def test_simple(helpers, emitter_directory, maf_file, cellline_lookup_path):
    """ simple test """
    validate(helpers, emitter_directory, maf_file, cellline_lookup_path)
