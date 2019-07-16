
import os
import contextlib
import pytest
import shutil
from transform.tcga.expression import transform, make_parallel_workstream


@pytest.fixture
def source_path(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/expression/transcript-level/TCGA_ACC_tpm.tsv.gz')


@pytest.fixture
def id_map_file(request):
    """ get the full path of the test id map file """
    return os.path.join(request.fspath.dirname, 'source/tcga/expression/transcript-level/TCGA_ID_MAP.csv')


@pytest.fixture
def source_wildcard(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/expression/transcript-level/*_tpm.tsv.gz')


def validate(helpers, source_path, id_map_file, emitter_directory):
    """ run xform and test results"""
    gene_expression_file = os.path.join(emitter_directory, 'ACC.GeneExpression.Vertex.json.gz')
    transcript_expression_file = os.path.join(emitter_directory, 'ACC.TranscriptExpression.Vertex.json.gz')

    agg_edge_file = os.path.join(emitter_directory, 'ACC.Aliquot_GeneExpressions_GeneExpression.Edge.json.gz')
    gaa_edge_file = os.path.join(emitter_directory, 'ACC.GeneExpression_Aliquot_Aliquot.Edge.json.gz')
    att_edge_file = os.path.join(emitter_directory, 'ACC.Aliquot_TranscriptExpressions_TranscriptExpression.Edge.json.gz')
    taa_edge_file = os.path.join(emitter_directory, 'ACC.TranscriptExpression_Aliquot_Aliquot.Edge.json.gz')

    all_files = [gene_expression_file, transcript_expression_file,
                 agg_edge_file, gaa_edge_file, att_edge_file, taa_edge_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(source_path=source_path, id_map_file=id_map_file, emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])


def test_simple(helpers, source_path, id_map_file, emitter_directory):
    """ just run validate"""
    validate(helpers, source_path, id_map_file, emitter_directory)


def test_make_parallel_workstream(source_wildcard):
    """ ensure main utility works """
    cmd = make_parallel_workstream(source_path=source_wildcard, jobs=5, dry_run=True)
    assert '--jobs 5' in cmd
    assert os.path.isfile('/tmp/tcga_expression_transform.txt'), 'should have created tcga_expression_transform.txt'
