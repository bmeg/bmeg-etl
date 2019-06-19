import pytest
import os
import contextlib
import shutil

from transform.go.go_obo2schema import transform


@pytest.fixture
def obo_file(request):
    return os.path.join(request.fspath.dirname, 'source/go/go.obo')


def test_simple(helpers, emitter_directory, obo_file):
    gene_ontology_term_file = os.path.join(emitter_directory, "GeneOntologyTerm.Vertex.json.gz")
    parents_edge_file = os.path.join(emitter_directory, "parent_terms.Edge.json.gz")
    children_edge_file = os.path.join(emitter_directory, "child_terms.Edge.json.gz")

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # run transform
    transform(obo_file=obo_file,
              emitter_directory=emitter_directory)

    # ratify
    count = helpers.assert_vertex_file_valid(gene_ontology_term_file)
    assert count == 6, "unexpected number of gene ontology terms emitted"
    helpers.assert_edge_file_valid(parents_edge_file)
    helpers.assert_edge_file_valid(children_edge_file)
    helpers.assert_edge_joins_valid([gene_ontology_term_file, parents_edge_file, children_edge_file])
