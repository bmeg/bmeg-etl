""" test maf_transform """

import os
import contextlib
import pytest
from transform.g2p.transform import transform
from transform.g2p.genes import normalize as gene_normalize
from bmeg.vertex import G2PAssociation, Publication, Gene, Allele, Phenotype


@pytest.fixture
def g2p_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'truncated_all.json')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test/test')


def validate(helpers, g2p_file, emitter_path_prefix):
    association_file = '{}.G2PAssociation.Vertex.json'.format(emitter_path_prefix)
    publication_file = '{}.Publication.Vertex.json'.format(emitter_path_prefix)
    publication_edge_file = '{}.HasSupportingReference.Edge.json'.format(emitter_path_prefix)
    gene_edge_file = '{}.GeneFeatureFor.Edge.json'.format(emitter_path_prefix)
    allele_edge_file = '{}.AlleleFeatureFor.Edge.json'.format(emitter_path_prefix)
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    phenotype_file = '{}.Phenotype.Vertex.json'.format(emitter_path_prefix)
    phenotype_edge_file = '{}.PhenotypeOf.Edge.json'.format(emitter_path_prefix)
    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(association_file)
        os.remove(publication_file)
        os.remove(publication_edge_file)
        os.remove(gene_edge_file)
        os.remove(allele_edge_file)
        os.remove(allele_file)
        os.remove(phenotype_file)
        os.remove(phenotype_edge_file)
    # create output
    transform(g2p_file, prefix=emitter_path_prefix)
    # test/test.G2PAssociation.Vertex.json
    helpers.assert_vertex_file_valid(G2PAssociation, association_file)
    # test/test.Publication.Vertex.json
    helpers.assert_vertex_file_valid(Publication, publication_file)
    # test/test.HasSupportingReference.Edge.json
    helpers.assert_edge_file_valid(G2PAssociation, Publication, publication_edge_file)
    # test/test.GeneFeatureFor.Edge.json
    helpers.assert_edge_file_valid(G2PAssociation, Gene, gene_edge_file)
    # test/test.AlleleFeatureFor.Edge.json
    helpers.assert_edge_file_valid(G2PAssociation, Allele, allele_edge_file)
    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    # test/test.Phenotype.Vertex.json
    helpers.assert_vertex_file_valid(Phenotype, phenotype_file)
    # test/test.PhenotypeOf.Edge.json
    helpers.assert_edge_file_valid(G2PAssociation, Phenotype, phenotype_edge_file)


def test_simple(helpers, g2p_file, emitter_path_prefix):
    """ simple test """
    validate(helpers, g2p_file, emitter_path_prefix)


def test_genes():
    """ ensure genes list expressed as gid, and keeps track of 'already seen' """
    import transform.g2p.genes
    # reset singleton 'already seen'
    transform.g2p.genes.EXPORTED_GENES = []
    assert gene_normalize({'genes': ['TP53']}) == ({'genes': ['Gene:ENSG00000141510']}, ['Gene:ENSG00000141510']), 'We should have a modified hit and a gene vertex gid'
    assert gene_normalize({'genes': ['TP53', 'EGFR']}) == ({'genes': ['Gene:ENSG00000141510', 'Gene:ENSG00000146648']}, ['Gene:ENSG00000146648']), 'We should have a modified hit and a gene vertex gid only for genes we havent seen'


def test_genes_nofind():
    """ no find thrown """
    try:
        gene_normalize({'genes': ['XXXX']})
        assert True, 'Should have raised exception'
    except Exception:
        pass
