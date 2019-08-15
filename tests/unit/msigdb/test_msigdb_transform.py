import contextlib
import os
import shutil
import pytest
from transform.msigdb.transform import transform


@pytest.fixture
def input_xml(request):
    return os.path.join(request.fspath.dirname, 'source/msigdb/msigdb_v6.2.xml')


def validate(helpers, emitter_directory, input_xml):
    """ run xform and test results"""
    geneset_file = os.path.join(emitter_directory, 'GeneSet.Vertex.json.gz')

    gs_gene_edge = os.path.join(emitter_directory, 'GeneSet_Genes_Gene.Edge.json.gz')
    gene_gs_edge = os.path.join(emitter_directory, 'Gene_GeneSets_GeneSet.Edge.json.gz')
    gs_publication_edge = os.path.join(emitter_directory, 'GeneSet_Publications_Publication.Edge.json.gz')
    publication_gs_edge = os.path.join(emitter_directory, 'Publication_GeneSets_GeneSet.Edge.json.gz')

    all_files = [
        # vertices
        geneset_file,
        # edges
        gs_gene_edge, gene_gs_edge,
        gs_publication_edge, publication_gs_edge
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(
        input_path=input_xml,
        emitter_prefix=None,
        emitter_directory=emitter_directory
    )

    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=["Gene", "Publication"]
    )


def test_simple(helpers, emitter_directory, input_xml):

    validate(helpers, emitter_directory, input_xml)
