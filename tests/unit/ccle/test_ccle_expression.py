import contextlib
import os
import shutil
import pytest
from transform.ccle.ccle_expression import transform_tpm, transform_gene_tpm


@pytest.fixture
def project_lookup_path(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/cellline_project_lookup.tsv')


@pytest.fixture
def gene_tpm_file(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_depMap_19Q1_TPM.csv')


@pytest.fixture
def tpm_file(request):
    return os.path.join(request.fspath.dirname, 'source/ccle/CCLE_depMap_19Q1_TPM_transcripts.csv')


def validate(helpers, emitter_directory, gene_tpm_file, tpm_file, project_lookup_path):
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
        shutil.rmtree(emitter_directory)

    # create output
    transform_gene_tpm(path=gene_tpm_file,
                       project_lookup_path=project_lookup_path,
                       emitter_directory=emitter_directory)

    transform_tpm(path=tpm_file,
                  project_lookup_path=project_lookup_path,
                  emitter_directory=emitter_directory)

    # ratify
    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(all_files, exclude_labels=['Aliquot'])


def test_simple(helpers, emitter_directory, gene_tpm_file, tpm_file, project_lookup_path):
    """ just run validate"""
    validate(helpers, emitter_directory, gene_tpm_file, tpm_file, project_lookup_path)
