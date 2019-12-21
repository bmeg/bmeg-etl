import os
import contextlib
import pytest
import shutil
import json

from transform.dgidb.transform import transform
from bmeg.ioutils import reader


@pytest.fixture
def interactions_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'source/dgidb/interactions.tsv')


def validate(helpers, interactions_file, emitter_directory):
    association_file = os.path.join(emitter_directory, 'G2PAssociation.Vertex.json.gz')
    compound_file = os.path.join(emitter_directory, 'Compound.Vertex.json.gz')

    pub_edge_file = os.path.join(emitter_directory, 'G2PAssociation_Publications_Publication.Edge.json.gz')
    pub_g2p_edge_file = os.path.join(emitter_directory, 'Publication_G2PAssociations_G2PAssociation.Edge.json.gz')

    gene_edge_file = os.path.join(emitter_directory, 'G2PAssociation_Genes_Gene.Edge.json.gz')
    gene_g2p_edge_file = os.path.join(emitter_directory, 'Gene_G2PAssociations_G2PAssociation.Edge.json.gz')

    compound_edge_file = os.path.join(emitter_directory, 'G2PAssociation_Compounds_Compound.Edge.json.gz')
    compound_g2p_edge_file = os.path.join(emitter_directory, 'Compound_G2PAssociations_G2PAssociation.Edge.json.gz')

    all_files = [
        association_file, compound_file,
        pub_edge_file, pub_g2p_edge_file, gene_edge_file, gene_g2p_edge_file,
        compound_edge_file, compound_g2p_edge_file
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(interactions_file,
              emitter_directory=emitter_directory)

    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
        exclude_labels=['Publication', 'Gene']
    )

    count = 0
    with reader(pub_edge_file) as f:
        for line in f:
            e = json.loads(line)
            assert e['to'] != 'Publication:ncbi.nlm.nih.gov/pubmed/'
            count += 1
    assert count == 16


def test_simple(helpers, interactions_file, emitter_directory):
    """ simple test """
    validate(helpers, interactions_file, emitter_directory)
