""" test maf_transform """

import os
import contextlib
import pytest
import shutil
from transform.g2p.transform import transform
from transform.g2p.genes import normalize as gene_normalize


@pytest.fixture
def g2p_file(request):
    """ get the full path of the test fixture """
    return os.path.join(request.fspath.dirname, 'truncated_all.json')


def validate(helpers, g2p_file, emitter_directory):
    association_file = os.path.join(emitter_directory, 'G2pAssociation.Vertex.json.gz')
    publication_edge_file = os.path.join(emitter_directory, 'publications.Edge.json.gz')
    gene_edge_file = os.path.join(emitter_directory, 'genes.Edge.json.gz')
    allele_edge_file = os.path.join(emitter_directory, 'alleles.Edge.json.gz')
    allele_file = os.path.join(emitter_directory, 'Allele.Vertex.json.gz')
    phenotype_file = os.path.join(emitter_directory, 'Phenotype.Vertex.json.gz')
    phenotype_edge_file = os.path.join(emitter_directory, 'phenotypes.Edge.json.gz')
    deadletter_file = os.path.join(emitter_directory, 'Deadletter.Vertex.json.gz')
    genomic_feature_file = os.path.join(emitter_directory, 'GenomicFeature.Vertex.json.gz')
    genomic_feature_edge_file = os.path.join(emitter_directory, 'genomic_features.Edge.json.gz')
    allele_in_edge_file = os.path.join(emitter_directory, 'gene.Edge.json.gz')
    environment_in_edge_file = os.path.join(emitter_directory, 'compounds.Edge.json.gz')
    environment_file = os.path.join(emitter_directory, 'Compound.Vertex.json.gz')
    g2p_edges_file = os.path.join(emitter_directory, 'g2p_associations.Edge.json.gz')

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(g2p_file, prefix=emitter_directory)

    # test/G2PAssociation.Vertex.json
    association_count = helpers.assert_vertex_file_valid(association_file)
    assert association_count == 29, 'There should be 29 associations'
    # test/publications.Edge.json
    helpers.assert_edge_file_valid(publication_edge_file)
    # genes.Edge.json.gz
    gene_count = helpers.assert_edge_file_valid(gene_edge_file)
    assert 52 == gene_count, 'There should be 52 edges to genes'

    from_tos = {}
    to_froms = {}
    # gene_edge_file = 'outputs/g2p/genes.Edge.json.gz'
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

    # test/alleles.Edge.json
    allele_count = helpers.assert_edge_file_valid(allele_edge_file)
    assert 78 == allele_count, 'There should be 78 alleles {}'.format(allele_edge_file)
    # test/Allele.Vertex.json
    helpers.assert_vertex_file_valid(allele_file)
    # test/Phenotype.Vertex.json
    phenotype_count = helpers.assert_vertex_file_valid(phenotype_file)
    assert phenotype_count == 209, 'There should be 209 phenotypes'
    # test/phenotypes.Edge.json
    has_phenotype_count = helpers.assert_edge_file_valid(phenotype_edge_file)
    assert has_phenotype_count == 295, 'There should be 295 has_phenotype edges'
    # test/Deadletter.Vertex.json
    helpers.assert_vertex_file_valid(deadletter_file)
    # test/GenomicFeature.Vertex.json
    helpers.assert_vertex_file_valid(genomic_feature_file)
    # test/Environment.Vertex.json
    compound_count = helpers.assert_vertex_file_valid(environment_file)
    assert compound_count == 46, 'There should be 46 compounds'
    # test/GenomicFeature.Vertex.json
    genomic_feature_edge_count = helpers.assert_edge_file_valid(genomic_feature_edge_file)
    assert genomic_feature_edge_count == 25, 'There should be 25 edges between G2PAssociation and GenomicFeature'
    # test/gene.Edge.json
    helpers.assert_edge_file_valid(allele_in_edge_file)
    # test/compounds.Edge.json
    helpers.assert_edge_file_valid(environment_in_edge_file)
    # test/g2p_associations.Edge.json
    helpers.assert_edge_file_valid(g2p_edges_file)

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        [association_file, publication_edge_file, gene_edge_file, allele_edge_file, allele_file,
         phenotype_file, phenotype_edge_file, genomic_feature_file, genomic_feature_edge_file,
         allele_in_edge_file, environment_in_edge_file, environment_file, g2p_edges_file],
        exclude_labels=['Publication', 'Gene']
    )


def test_simple(helpers, g2p_file, emitter_directory):
    """ simple test """
    validate(helpers, g2p_file, emitter_directory)


def test_genes():
    """ ensure genes list expressed as gid, and keeps track of 'already seen' """
    import transform.g2p.genes
    # reset singleton 'already seen'
    transform.g2p.genes.EXPORTED_GENES = []
    normalized = gene_normalize({'genes': ['TP53']})
    assert normalized == ({'genes': ['ENSG00000141510']}, ['ENSG00000141510'], []), 'We should have a modified hit and a gene vertex gid'
    normalized = gene_normalize({'genes': ['TP53', 'EGFR']})
    print(normalized)
    assert normalized[0] == {'genes': ['ENSG00000146648']}, 'should return both genes'
    assert 'ENSG00000146648' in normalized[1] and 'ENSG00000141510' not in normalized[1], 'should return both genes'


def test_genes_nofind():
    """ no find thrown """
    try:
        gene_normalize({'genes': ['XXXX']})
        assert True, 'Should have raised exception'
    except Exception:
        pass
