import contextlib
import os
import shutil
import pytest

from transform.gtex.expression import transform


@pytest.fixture
def gct_file(request):
    return os.path.join(request.fspath.dirname, 'source/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz')


@pytest.fixture
def project_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/gtex/project_lookup.tsv')


def validate(helpers, emitter_directory, gct_file, project_lookup_path):
    """ run xform and test results"""
    gene_expression_file = os.path.join(emitter_directory, 'GeneExpression.Vertex.json.gz')

    gene_expressions_edge_file = os.path.join(emitter_directory, 'gene_expressions.Edge.json.gz')
    aliquot_edge_file = os.path.join(emitter_directory, 'aliquot.Edge.json.gz')

    all_files = [gene_expression_file,
                 gene_expressions_edge_file, aliquot_edge_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(gct_file=gct_file,
              project_lookup_path=project_lookup_path,
              emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])


def test_simple(helpers, emitter_directory, gct_file, project_lookup_path):
    """ just run validate"""
    validate(helpers, emitter_directory, gct_file, project_lookup_path)
