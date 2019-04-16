""" test maf_transform """

import os
import contextlib
import pytest
from transform.g2p.transform import transform
from transform.g2p.genes import normalize as gene_normalize
from bmeg.vertex import G2PAssociation, Publication, Gene, Allele, Phenotype, Deadletter, GenomicFeature, Compound


@pytest.fixture
def g2p_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'truncated_all.json')


@pytest.fixture
def emitter_path_prefix(request):
    """ get the full path of the test output """
    return os.path.join(request.fspath.dirname, 'test')


def validate(helpers, g2p_file, emitter_path_prefix):
    association_file = os.path.join(emitter_path_prefix, 'G2PAssociation.Vertex.json.gz')
    publication_edge_file = os.path.join(emitter_path_prefix, 'HasSupportingReference.Edge.json.gz')
    gene_edge_file = os.path.join(emitter_path_prefix, 'HasGeneFeature.Edge.json.gz')
    allele_edge_file = os.path.join(emitter_path_prefix, 'HasAlleleFeature.Edge.json.gz')
    allele_file = os.path.join(emitter_path_prefix, 'Allele.Vertex.json.gz')
    phenotype_file = os.path.join(emitter_path_prefix, 'Phenotype.Vertex.json.gz')
    phenotype_edge_file = os.path.join(emitter_path_prefix, 'HasPhenotype.Edge.json.gz')
    deadletter_file = os.path.join(emitter_path_prefix, 'Deadletter.Vertex.json.gz')
    genomic_feature_file = os.path.join(emitter_path_prefix, 'GenomicFeature.Vertex.json.gz')
    genomic_feature_edge_file = os.path.join(emitter_path_prefix, 'HasGenomicFeature.Edge.json.gz')
    has_gene_edge_file = os.path.join(emitter_path_prefix, 'GenomicFeatureIn.Edge.json.gz')
    allele_in_edge_file = os.path.join(emitter_path_prefix, 'AlleleIn.Edge.json.gz')
    environment_in_edge_file = os.path.join(emitter_path_prefix, 'HasEnvironment.Edge.json.gz')
    environment_file = os.path.join(emitter_path_prefix, 'Compound.Vertex.json.gz')
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
        os.remove(genomic_feature_file)
        os.remove(genomic_feature_edge_file)
        os.remove(has_gene_edge_file)
        os.remove(allele_in_edge_file)
        os.remove(environment_in_edge_file)
        os.remove(environment_file)

    # create output
    transform(g2p_file, prefix=emitter_path_prefix)
    # test/test.G2PAssociation.Vertex.json
    association_count = helpers.assert_vertex_file_valid(G2PAssociation, association_file)
    assert association_count == 29, 'There should be 29 associations'
    # test/test.HasSupportingReference.Edge.json
    helpers.assert_edge_file_valid(G2PAssociation, Publication, publication_edge_file)
    # HasGeneFeature.Edge.json.gz
    gene_count = helpers.assert_edge_file_valid(G2PAssociation, Gene, gene_edge_file)
    assert 46 == gene_count, 'There should be 46 genes to association'

    from_tos = {}
    to_froms = {}
    # gene_edge_file = 'outputs/g2p/HasGeneFeature.Edge.json.gz'
    for edge in helpers.parse_edge_file(gene_edge_file):
        from_ = edge['from']
        to_ = edge['to']

        tos = from_tos.get(from_, set())
        tos.add(to_)
        from_tos[from_] = tos

        froms = to_froms.get(to_, set())
        froms.add(from_)
        to_froms[to_] = froms

    multiple_count = 0
    for from_ in from_tos:
        if len(from_tos[from_]) > 1:
            multiple_count += 1
    assert multiple_count > 0, 'at least one association should have multiple genes {}'.format(from_tos)

    multiple_count = 0
    for to_ in to_froms:
        if len(to_froms[to_]) > 1:
            multiple_count += 1
    assert multiple_count > 0, 'at least one gene should have multiple associations {}'.format(to_froms)

    # test/test.AlleleFeatureFor.Edge.json
    allele_count = helpers.assert_edge_file_valid(G2PAssociation, Allele, allele_edge_file)
    assert 39 == allele_count, 'There should be 39 alleles {}'.format(allele_edge_file)
    # test/test.Allele.Vertex.json
    helpers.assert_vertex_file_valid(Allele, allele_file)
    # test/test.Phenotype.Vertex.json
    phenotype_count = helpers.assert_vertex_file_valid(Phenotype, phenotype_file)
    assert phenotype_count == 209, 'There should be 209 phenotypes'
    # test/test.HasPhenotype.Edge.json
    has_phenotype_count = helpers.assert_edge_file_valid(G2PAssociation, Phenotype, phenotype_edge_file)
    assert has_phenotype_count == 295, 'There should be 295 has_phenotype edges'
    # test/test.Deadletter.Vertex.json
    helpers.assert_vertex_file_valid(Deadletter, deadletter_file)
    # test/test.GenomicFeature.Vertex.json
    helpers.assert_vertex_file_valid(GenomicFeature, genomic_feature_file)
    # test/test.Environment.Vertex.json
    compound_count = helpers.assert_vertex_file_valid(Compound, environment_file)
    assert compound_count == 46, 'There should be 46 compounds'
    # test/test.GenomicFeature.Vertex.json
    genomic_feature_edge_count = helpers.assert_edge_file_valid(G2PAssociation, GenomicFeature, genomic_feature_edge_file)
    assert genomic_feature_edge_count == 14, 'There should be 14 edges between G2PAssociation and GenomicFeature'
    # test/test.HasGene.Edge.json
    helpers.assert_edge_file_valid(GenomicFeature, Gene, has_gene_edge_file)
    # test/test.AlleleIn.Edge.json
    helpers.assert_edge_file_valid(Allele, Gene, allele_in_edge_file)
    # test/test.HasEnvironment.Edge.json
    helpers.assert_edge_file_valid(G2PAssociation, Compound, environment_in_edge_file)
    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        [association_file, publication_edge_file, gene_edge_file, allele_edge_file, allele_file, phenotype_file, phenotype_edge_file, deadletter_file, genomic_feature_file, genomic_feature_edge_file, has_gene_edge_file, allele_in_edge_file, environment_in_edge_file, environment_file],
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
    assert gene_normalize({'genes': ['TP53']}) == ({'genes': {'ENSG00000141510'}}, ['ENSG00000141510'], []), 'We should have a modified hit and a gene vertex gid'
    normalized = gene_normalize({'genes': ['TP53', 'EGFR']})
    assert normalized[0] == {'genes': {'ENSG00000141510', 'ENSG00000146648'}}, 'should return both genes'
    assert 'ENSG00000146648' in normalized[1] and 'ENSG00000141510' in normalized[1], 'should return both genes'


def test_genes_nofind():
    """ no find thrown """
    try:
        gene_normalize({'genes': ['XXXX']})
        assert True, 'Should have raised exception'
    except Exception:
        pass
