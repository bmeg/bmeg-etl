import os
import contextlib
import pytest
import shutil
from transform.tcga.gistic2_cna import transform, make_parallel_workstream


@pytest.fixture
def source_path(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/gistic2-firehose/TCGA-ACC_all_thresholded.by_genes.txt')


@pytest.fixture
def source_wildcard(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/gistic2-firehose/*_all_thresholded.by_genes.txt')


@pytest.fixture
def id_lookup_path(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdc/id_lookup.tsv')


@pytest.fixture
def project_lookup_path(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/gdc/project_lookup.tsv')


def validate(helpers, emitter_directory, source_path, id_lookup_path, project_lookup_path):
    """ run xform and test results"""
    cna_file = os.path.join(emitter_directory, 'TCGA-ACC.CopyNumberAlteration.Vertex.json.gz')
    cna_edge_file = os.path.join(emitter_directory, 'TCGA-ACC.copy_number_alterations.Edge.json.gz')
    aliquot_edge_file = os.path.join(emitter_directory, 'TCGA-ACC.aliquot.Edge.json.gz')

    all_files = [cna_file, cna_edge_file, aliquot_edge_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(source_path=source_path,
              id_lookup_path=id_lookup_path,
              project_lookup_path=project_lookup_path,
              emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])


def test_simple(helpers, emitter_directory, source_path, id_lookup_path, project_lookup_path):
    """ just run validate"""
    validate(helpers, emitter_directory, source_path, id_lookup_path, project_lookup_path)


def test_make_parallel_workstream(source_wildcard):
    """ ensure main utility works """
    cmd = make_parallel_workstream(source_path=source_wildcard, jobs=5, dry_run=True)
    assert '--jobs 5' in cmd
    assert os.path.isfile('/tmp/tcga_gistic2_transform.txt'), 'should have created tcga_gistic2_transform.txt'
