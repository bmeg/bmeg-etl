
import os
import contextlib
import pytest
from transform.ccle.ccle_expression import transform_tpm, transform_gene_tpm


@pytest.fixture
def gene_tpm_file(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_depMap_19Q1_TPM.csv')


@pytest.fixture
def tpm_file(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_depMap_19Q1_TPM_transcripts.csv')


def validate(helpers, gene_tpm_file, tpm_file, emitter_directory):
    """ run xform and test results"""
    gene_expression_file = os.path.join(emitter_directory, 'GeneExpression.Vertex.json.gz')
    transcript_expression_file = os.path.join(emitter_directory, 'TranscriptExpression.Vertex.json.gz')

    gene_expressions_edge_file = os.path.join(emitter_directory, 'gene_expressions.Edge.json.gz')
    transcript_expressions_edge_file = os.path.join(emitter_directory, 'transcript_expressions.Edge.json.gz')
    aliquot_edge_file = os.path.join(emitter_directory, 'aliquot.Edge.json.gz')

    all_files = [gene_expression_file, transcript_expression_file,
                 gene_expressions_edge_file, transcript_expressions_edge_file,
                 aliquot_edge_file]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        for f in all_files:
            os.remove(f)

    # create output
    transform_gene_tpm(path=gene_tpm_file, emitter_directory=emitter_directory)
    transform_tpm(path=tpm_file, emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])


def test_simple(helpers, gene_tpm_file, tpm_file, emitter_directory):
    """ just run validate"""
    validate(helpers, gene_tpm_file, tpm_file, emitter_directory)
