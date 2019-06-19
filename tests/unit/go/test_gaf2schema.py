import pytest
import os
import contextlib
import shutil

from transform.go.go_gaf2schema import transform


@pytest.fixture
def gaf_file(request):
    return os.path.join(request.fspath.dirname, 'source/go/goa_human.gaf.gz')


@pytest.fixture
def id_map_file(request):
    return os.path.join(request.fspath.dirname, 'source/go/HUMAN_9606_idmapping.dat.gz')


def test_simple(helpers, emitter_directory, gaf_file, id_map_file):
    genes_edge_file = os.path.join(emitter_directory, "genes.Edge.json.gz")
    gene_ontology_terms_edge_file = os.path.join(emitter_directory, "gene_ontology_terms.Edge.json.gz")

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # run transform
    transform(gaf_file=gaf_file,
              id_map_file=id_map_file,
              emitter_directory=emitter_directory)

    # ratify
    helpers.assert_edge_file_valid(genes_edge_file)
    helpers.assert_edge_file_valid(gene_ontology_terms_edge_file)
