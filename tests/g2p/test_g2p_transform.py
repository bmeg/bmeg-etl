""" test maf_transform """

import os
import contextlib
import pytest
from transform.g2p.transform import transform
from transform.g2p.genes import normalize as gene_normalize
from bmeg.vertex import G2PAssociation, Publication, Gene, Allele, Phenotype, Deadletter, MinimalAllele


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
    publication_edge_file = '{}.HasSupportingReference.Edge.json'.format(emitter_path_prefix)
    gene_edge_file = '{}.HasGeneFeature.Edge.json'.format(emitter_path_prefix)
    allele_edge_file = '{}.HasAlleleFeature.Edge.json'.format(emitter_path_prefix)
    allele_file = '{}.Allele.Vertex.json'.format(emitter_path_prefix)
    phenotype_file = '{}.Phenotype.Vertex.json'.format(emitter_path_prefix)
    phenotype_edge_file = '{}.HasPhenotype.Edge.json'.format(emitter_path_prefix)
    deadletter_file = '{}.Deadletter.Vertex.json'.format(emitter_path_prefix)
    minimal_allele_file = '{}.MinimalAllele.Vertex.json'.format(emitter_path_prefix)
    minimal_allele_edge_file = '{}.HasMinimalAlleleFeature.Edge.json'.format(emitter_path_prefix)
    has_gene_edge_file = '{}.HasGene.Edge.json'.format(emitter_path_prefix)
    allele_in_edge_file = '{}.AlleleIn.Edge.json'.format(emitter_path_prefix)
    # remove output
    with contextlib.suppress(FileNotFoundError):
        os.remove(association_file)
        os.remove(publication_edge_file)
        os.remove(gene_edge_file)
        os.remove(allele_edge_file)
        os.remove(allele_file)
        os.remove(phenotype_file)
        os.remove(phenotype_edge_file)
        os.remove(deadletter_file)
        os.remove(minimal_allele_file)
        os.remove(minimal_allele_edge_file)
        os.remove(has_gene_edge_file)
        os.remove(allele_in_edge_file)

    # create output
    transform(g2p_file, prefix=emitter_path_prefix)
    # test/test.G2PAssociation.Vertex.json
    helpers.assert_vertex_file_valid(G2PAssociation, association_file)
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
    # test/test.Deadletter.Vertex.json
    helpers.assert_vertex_file_valid(Deadletter, deadletter_file)
    # test/test.MinimalAllele.Vertex.json
    helpers.assert_vertex_file_valid(MinimalAllele, minimal_allele_file)
    # test/test.MinimalAllele.Vertex.json
    helpers.assert_edge_file_valid(G2PAssociation, MinimalAllele, minimal_allele_edge_file)
    # test/test.HasGene.Edge.json
    helpers.assert_edge_file_valid(MinimalAllele, Gene, has_gene_edge_file)
    # test/test.AlleleIn.Edge.json
    helpers.assert_edge_file_valid(Allele, Gene, allele_in_edge_file)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        [association_file, publication_edge_file, gene_edge_file, allele_edge_file, allele_file, phenotype_file, phenotype_edge_file, deadletter_file, minimal_allele_file, minimal_allele_edge_file, has_gene_edge_file, allele_in_edge_file],
        exclude_labels=['Publication', 'Gene']
    )


def test_simple(helpers, g2p_file, emitter_path_prefix):
    """ simple test """
    validate(helpers, g2p_file, emitter_path_prefix)


def test_genes():
    """ ensure genes list expressed as gid, and keeps track of 'already seen' """
    import transform.g2p.genes
    # reset singleton 'already seen'
    transform.g2p.genes.EXPORTED_GENES = []
    assert gene_normalize({'genes': ['TP53']}) == ({'genes': ['Gene:ENSG00000141510']}, ['Gene:ENSG00000141510'], []), 'We should have a modified hit and a gene vertex gid'
    assert gene_normalize({'genes': ['TP53', 'EGFR']}) == ({'genes': ['Gene:ENSG00000141510', 'Gene:ENSG00000146648']}, ['Gene:ENSG00000146648'], []), 'We should have a modified hit and a gene vertex gid only for genes we havent seen'


def test_genes_nofind():
    """ no find thrown """
    try:
        gene_normalize({'genes': ['XXXX']})
        assert True, 'Should have raised exception'
    except Exception:
        pass
