
import os
import contextlib
import pytest
from transform.tcga.gistic2_cna import transform, make_parallel_workstream
from bmeg.vertex import CopyNumberAlteration, Aliquot


@pytest.fixture
def source_path(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/gistic2-firehose/TCGA-ACC_all_thresholded.by_genes.txt')


@pytest.fixture
def source_wildcard(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/gistic2-firehose/*_all_thresholded.by_genes.txt')


def validate(helpers, source_path, emitter_directory):
    """ run xform and test results"""
    cna_file = os.path.join(emitter_directory, 'TCGA-ACC.CopyNumberAlteration.Vertex.json')
    cna_of_file = os.path.join(emitter_directory, 'TCGA-ACC.CopyNumberAlterationOf.Edge.json')

    all_files = [cna_file, cna_of_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform(source_path=source_path, emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(CopyNumberAlteration, cna_file)
    helpers.assert_edge_file_valid(CopyNumberAlteration, Aliquot, cna_of_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])


def test_simple(helpers, source_path, emitter_directory):
    """ just run validate"""
    validate(helpers, source_path, emitter_directory)


def test_make_parallel_workstream(source_wildcard):
    """ ensure main utility works """
    cmd = make_parallel_workstream(source_path=source_wildcard, jobs=5, dry_run=True)
    assert '--jobs 5' in cmd
    assert os.path.isfile('/tmp/tcga_gistic2_transform.txt'), 'should have created tcga_gistic2_transform.txt'
