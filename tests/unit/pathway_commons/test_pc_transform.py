import contextlib
import os
import shutil
import pytest
from transform.pathway_commons.transform import transform


@pytest.fixture
def sif(request):
    return os.path.join(request.fspath.dirname, 'source/pathway_commons/R-HSA-111448.sif')


@pytest.fixture
def pathways(request):
    return os.path.join(request.fspath.dirname, 'source/pathway_commons/pathways.R-HSA-111448.txt')


def validate(helpers, emitter_directory, sif, pathways):
    """ run xform and test results"""
    interaction_file = os.path.join(emitter_directory, 'Interaction.Vertex.json.gz')
    pathway_file = os.path.join(emitter_directory, 'Pathway.Vertex.json.gz')

    pathway_gene_edge = os.path.join(emitter_directory, 'Pathway_Genes_Gene.Edge.json.gz')
    gene_pathway_edge = os.path.join(emitter_directory, 'Gene_Pathways_Pathway.Edge.json.gz')
    pathway_interaction_edge = os.path.join(emitter_directory, 'Pathway_Interactions_Interaction.Edge.json.gz')
    interaction_pathway_edge = os.path.join(emitter_directory, 'Interaction_Pathways_Pathway.Edge.json.gz')
    interaction_gene_edge = os.path.join(emitter_directory, 'Interaction_Genes_Gene.Edge.json.gz')
    gene_interaction_edge = os.path.join(emitter_directory, 'Gene_Interactions_Interaction.Edge.json.gz')
    interaction_publication_edge = os.path.join(emitter_directory, 'Interaction_Publications_Publication.Edge.json.gz')
    publication_interaction_edge = os.path.join(emitter_directory, 'Publication_Interactions_Interaction.Edge.json.gz')

    all_files = [
        # vertices
        interaction_file, pathway_file,
        # edges
        pathway_gene_edge, gene_pathway_edge,
        pathway_interaction_edge, interaction_pathway_edge,
        interaction_gene_edge, gene_interaction_edge,
        interaction_publication_edge, publication_interaction_edge
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(
        sif_file=sif,
        pathways_file=pathways,
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


def test_simple(helpers, emitter_directory, sif, pathways):

    validate(helpers, emitter_directory, sif, pathways)
