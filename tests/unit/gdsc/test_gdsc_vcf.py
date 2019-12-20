import pytest
from transform.gdsc.vcf_transform import transform
from bmeg.ioutils import reader

import os
import contextlib
import shutil
import json


@pytest.fixture
def vcf_dir(request):
    return os.path.join(request.fspath.dirname, 'source/gdsc/vcfs/*')


@pytest.fixture
def cellline_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_id_lookup.tsv')


def validate(helpers, emitter_directory, vcf_dir, cellline_lookup_path):
    allele_file = os.path.join(emitter_directory, 'Allele.Vertex.json.gz')
    callset_file = os.path.join(emitter_directory, 'SomaticCallset.Vertex.json.gz')

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
    transform(
        vcf_dir=vcf_dir,
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
            assert "GDSC" in callset['gid'], 'gid should contain GDSC'
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


def test_simple(helpers, emitter_directory, vcf_dir, cellline_lookup_path):
    """ simple test """
    validate(helpers, emitter_directory, vcf_dir, cellline_lookup_path)
