
import os
import contextlib
import pytest
from transform.tcga.expression import transform, make_parallel_workstream
from bmeg.vertex import Expression, Aliquot


@pytest.fixture
def source_path(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/expression/transcript-level/TCGA_ACC_tpm.tsv.gz')


@pytest.fixture
def id_map_file(request):
    """ get the full path of the test id map file """
    return os.path.join(request.fspath.dirname, 'source/tcga/expression/transcript-level/TCGA_ID_MAP.csv')


@pytest.fixture
def gene_map_file(request):
    """ get the full path of the test gene map file"""
    return os.path.join(request.fspath.dirname, 'source/ensembl/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.trans_gene.tsv')


@pytest.fixture
def source_wildcard(request):
    """ get the full path of the test input """
    return os.path.join(request.fspath.dirname, 'source/tcga/expression/transcript-level/*_tpm.tsv.gz')


def validate(helpers, source_path, id_map_file, gene_map_file, emitter_directory):
    """ run xform and test results"""
    expression_file = os.path.join(emitter_directory, 'ACC.Expression.Vertex.json.gz')
    expression_of_file = os.path.join(emitter_directory, 'ACC.ExpressionOf.Edge.json.gz')

    all_files = [expression_file, expression_of_file]
    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform(source_path=source_path, id_map_file=id_map_file, gene_map_file=gene_map_file, emitter_directory=emitter_directory)
    # ratify
    helpers.assert_vertex_file_valid(Expression, expression_file)
    helpers.assert_edge_file_valid(Expression, Aliquot, expression_of_file)
    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])


def test_simple(helpers, source_path, id_map_file, gene_map_file, emitter_directory):
    """ just run validate"""
    validate(helpers, source_path, id_map_file, gene_map_file, emitter_directory)


def test_make_parallel_workstream(source_wildcard):
    """ ensure main utility works """
    cmd = make_parallel_workstream(source_path=source_wildcard, jobs=5, dry_run=True)
    assert '--jobs 5' in cmd
    assert os.path.isfile('/tmp/tcga_expression_transform.txt'), 'should have created tcga_expression_transform.txt'
