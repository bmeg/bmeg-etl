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
    association_file = os.path.join(emitter_directory, 'G2PAssociation.Vertex.json.gz')
    phenotype_file = os.path.join(emitter_directory, 'Phenotype.Vertex.json.gz')
    deadletter_file = os.path.join(emitter_directory, 'Deadletter.Vertex.json.gz')
    genomic_feature_file = os.path.join(emitter_directory, 'GenomicFeature.Vertex.json.gz')
    allele_file = os.path.join(emitter_directory, 'Allele.Vertex.json.gz')
    compound_file = os.path.join(emitter_directory, 'Compound.Vertex.json.gz')

    pub_edge_file = os.path.join(emitter_directory, 'G2PAssociation_Publications_Publication.Edge.json.gz')
    pub_g2p_edge_file = os.path.join(emitter_directory, 'Publication_G2PAssociations_G2PAssociation.Edge.json.gz')

    gene_edge_file = os.path.join(emitter_directory, 'G2PAssociation_Genes_Gene.Edge.json.gz')
    gene_g2p_edge_file = os.path.join(emitter_directory, 'Gene_G2PAssociations_G2PAssociation.Edge.json.gz')

    allele_edge_file = os.path.join(emitter_directory, 'G2PAssociation_Alleles_Allele.Edge.json.gz')
    allele_g2p_edge_file = os.path.join(emitter_directory, 'Allele_G2PAssociations_G2PAssociation.Edge.json.gz')

    phenotype_edge_file = os.path.join(emitter_directory, 'G2PAssociation_Phenotypes_Phenotype.Edge.json.gz')
    phenotype_g2p_edge_file = os.path.join(emitter_directory, 'Phenotype_G2PAssociations_G2PAssociation.Edge.json.gz')

    genomic_feature_edge_file = os.path.join(emitter_directory, 'G2PAssociation_GenomicFeatures_GenomicFeature.Edge.json.gz')
    genomic_feature_g2p_edge_file = os.path.join(emitter_directory, 'GenomicFeature_G2PAssociations_G2PAssociation.Edge.json.gz')

    compound_edge_file = os.path.join(emitter_directory, 'G2PAssociation_Compounds_Compound.Edge.json.gz')
    compound_g2p_edge_file = os.path.join(emitter_directory, 'Compound_G2PAssociations_G2PAssociation.Edge.json.gz')

    allele_gene_edge_file = os.path.join(emitter_directory, 'Allele_Gene_Gene.Edge.json.gz')
    gene_allele_edge_file = os.path.join(emitter_directory, 'Gene_Alleles_Allele.Edge.json.gz')

    all_files = [
        association_file, phenotype_file, deadletter_file, genomic_feature_file, allele_file, compound_file,
        pub_edge_file, pub_g2p_edge_file, gene_edge_file, gene_g2p_edge_file, allele_edge_file, allele_g2p_edge_file,
        phenotype_edge_file, phenotype_g2p_edge_file, genomic_feature_edge_file, genomic_feature_g2p_edge_file,
        compound_edge_file, compound_g2p_edge_file, allele_gene_edge_file, gene_allele_edge_file
    ]

    # remove output
    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree(emitter_directory)

    # create output
    transform(g2p_file, prefix=emitter_directory)

    for f in all_files:
        if "Vertex.json.gz" in f:
            helpers.assert_vertex_file_valid(f)
        elif "Edge.json.gz" in f:
            helpers.assert_edge_file_valid(f)

    association_count = helpers.assert_vertex_file_valid(association_file)
    assert association_count == 29, 'There should be 29 associations'

    gene_count = helpers.assert_edge_file_valid(gene_edge_file)
    assert gene_count == 46, 'There should be 46 edges to genes {}'.format(gene_edge_file)
    gene_count = helpers.assert_edge_file_valid(gene_g2p_edge_file)
    assert gene_count == 46, 'There should be 46 edges to genes {}'.format(gene_g2p_edge_file)

    from_tos = {}
    to_froms = {}
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

    allele_count = helpers.assert_edge_file_valid(allele_edge_file)
    assert allele_count == 40, 'There should be 40 alleles {}'.format(allele_edge_file)
    allele_count = helpers.assert_edge_file_valid(allele_g2p_edge_file)
    assert allele_count == 40, 'There should be 40 alleles {}'.format(allele_g2p_edge_file)

    phenotype_count = helpers.assert_vertex_file_valid(phenotype_file)
    assert phenotype_count == 209, 'There should be 209 phenotypes'

    phenotype_edge_count = helpers.assert_edge_file_valid(phenotype_edge_file)
    assert phenotype_edge_count == 295, 'There should be 295 phenotype edges {}'.format(phenotype_edge_file)
    phenotype_edge_count = helpers.assert_edge_file_valid(phenotype_g2p_edge_file)
    assert phenotype_edge_count == 295, 'There should be 295 phenotype edges {}'.format(phenotype_g2p_edge_file)

    genomic_feature_edge_count = helpers.assert_edge_file_valid(genomic_feature_edge_file)
    assert genomic_feature_edge_count == 13, 'There should be 13 edges between G2PAssociation and GenomicFeature {}'.format(genomic_feature_edge_file)
    genomic_feature_edge_count = helpers.assert_edge_file_valid(genomic_feature_g2p_edge_file)
    assert genomic_feature_edge_count == 13, 'There should be 13 edges between G2PAssociation and GenomicFeature {}'.format(genomic_feature_g2p_edge_file)

    compound_count = helpers.assert_vertex_file_valid(compound_file)
    assert compound_count == 46, 'There should be 46 compounds'

    # validate vertex for all edges exist
    helpers.assert_edge_joins_valid(
        all_files,
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
    assert normalized[0] == {'genes': ['ENSG00000141510', 'ENSG00000146648']}, 'should return both genes'
    assert 'ENSG00000146648' in normalized[1] and 'ENSG00000141510' in normalized[1], 'should return both genes'


def test_genes_nofind():
    """ no find thrown """
    try:
        gene_normalize({'genes': ['XXXX']})
        assert True, 'Should have raised exception'
    except Exception:
        pass
