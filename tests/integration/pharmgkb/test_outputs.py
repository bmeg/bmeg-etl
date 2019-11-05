"""
test transform results in outputs/pharmgkb
"""
from bmeg.ioutils import reader
import os
import json
from collections import defaultdict

expected_files = """
    outputs/pharmgkb/Allele.Vertex.json.gz
    outputs/pharmgkb/Allele_G2PAssociations_G2PAssociation.Edge.json.gz
    outputs/pharmgkb/Allele_Gene_Gene.Edge.json.gz
    outputs/pharmgkb/Compound.Vertex.json.gz
    outputs/pharmgkb/Compound_G2PAssociations_G2PAssociation.Edge.json.gz
    outputs/pharmgkb/G2PAssociation.Vertex.json.gz
    outputs/pharmgkb/G2PAssociation_Alleles_Allele.Edge.json.gz
    outputs/pharmgkb/G2PAssociation_Compounds_Compound.Edge.json.gz
    outputs/pharmgkb/G2PAssociation_Genes_Gene.Edge.json.gz
    outputs/pharmgkb/G2PAssociation_GenomicFeatures_GenomicFeature.Edge.json.gz
    outputs/pharmgkb/G2PAssociation_Phenotypes_Phenotype.Edge.json.gz
    outputs/pharmgkb/G2PAssociation_Publications_Publication.Edge.json.gz
    outputs/pharmgkb/Gene_Alleles_Allele.Edge.json.gz
    outputs/pharmgkb/Gene_G2PAssociations_G2PAssociation.Edge.json.gz
    outputs/pharmgkb/GenomicFeature.Vertex.json.gz
    outputs/pharmgkb/GenomicFeature_G2PAssociations_G2PAssociation.Edge.json.gz
    outputs/pharmgkb/Phenotype.Vertex.json.gz
    outputs/pharmgkb/Phenotype_G2PAssociations_G2PAssociation.Edge.json.gz
    outputs/pharmgkb/Publication.Vertex.json.gz
    outputs/pharmgkb/Publication_G2PAssociations_G2PAssociation.Edge.json.gz
""".strip().split()

expected_file_counts = {
    'outputs/pharmgkb/Allele.Vertex.json.gz': 2911,
    'outputs/pharmgkb/Allele_G2PAssociations_G2PAssociation.Edge.json.gz': 18442,
    'outputs/pharmgkb/Allele_Gene_Gene.Edge.json.gz': 811,
    'outputs/pharmgkb/Compound.Vertex.json.gz': 551,
    'outputs/pharmgkb/Compound_G2PAssociations_G2PAssociation.Edge.json.gz': 5977,
    'outputs/pharmgkb/G2PAssociation.Vertex.json.gz': 4160,
    'outputs/pharmgkb/G2PAssociation_Alleles_Allele.Edge.json.gz': 18442,
    'outputs/pharmgkb/G2PAssociation_Compounds_Compound.Edge.json.gz': 5977,
    'outputs/pharmgkb/G2PAssociation_Genes_Gene.Edge.json.gz': 4332,
    'outputs/pharmgkb/G2PAssociation_GenomicFeatures_GenomicFeature.Edge.json.gz': 20428,
    'outputs/pharmgkb/G2PAssociation_Phenotypes_Phenotype.Edge.json.gz': 4148,
    'outputs/pharmgkb/G2PAssociation_Publications_Publication.Edge.json.gz': 10658,
    'outputs/pharmgkb/Gene_Alleles_Allele.Edge.json.gz': 811,
    'outputs/pharmgkb/Gene_G2PAssociations_G2PAssociation.Edge.json.gz': 4332,
    'outputs/pharmgkb/GenomicFeature.Vertex.json.gz': 3436,
    'outputs/pharmgkb/GenomicFeature_G2PAssociations_G2PAssociation.Edge.json.gz': 20428,
    'outputs/pharmgkb/Phenotype.Vertex.json.gz': 292,
    'outputs/pharmgkb/Phenotype_G2PAssociations_G2PAssociation.Edge.json.gz': 4148,
    'outputs/pharmgkb/Publication.Vertex.json.gz': 5105,
    'outputs/pharmgkb/Publication_G2PAssociations_G2PAssociation.Edge.json.gz': 10658,
}

normalized_files = """
    outputs/phenotype/normalized.G2PAssociation_Phenotypes_Phenotype.Edge.json.gz
    outputs/phenotype/normalized.Phenotype_G2PAssociations_G2PAssociation.Edge.json.gz
    outputs/compound/normalized.Compound_G2PAssociations_G2PAssociation.Edge.json.gz
    outputs/compound/normalized.G2PAssociation_Compounds_Compound.Edge.json.gz
""".strip().split()


def test_expected_files_exists():
    for f in expected_files:
        assert os.path.exists(f), f"{f} Should exist"


def test_expected_files_counts():
    for f in expected_file_counts:
        c = 0
        with reader(f) as ins:
            for l in ins:
                c += 1
        assert c == expected_file_counts[f], f"{f} expected:{expected_file_counts[f]} actual:{c}"


def test_clinical_ann_metadata_count():
    f = 'source/pharmgkb/clinical_ann_metadata.tsv'
    c = 0
    with reader(f) as ins:
        for l in ins:
            c += 1
    expected = expected_file_counts['outputs/pharmgkb/G2PAssociation.Vertex.json.gz']
    c = c - 1  # don't count header
    assert c == expected, f"clinical_ann_metadata actual: {c} expected: {expected}"


def test_clinical_ann_metadata_1447989723():
    f = 'outputs/pharmgkb/G2PAssociation.Vertex.json.gz'
    association_gid = None
    with reader(f) as ins:
        for l in ins:
            a = json.loads(l)  # assoc
            cam = json.loads(a['data']['source_document'])  # clinical_ann_metadata
            if cam['clinical_annotation_id'] == '1447989723':
                association_gid = a['gid']
                assert a['data']['evidence_label'] == '1A'
                assert a['data']['source'] == 'pharmgkb'
                assert a['data']['source_url'] == 'https://www.pharmgkb.org/clinicalAnnotation/1447989723'
    edges = defaultdict(list)
    for f in expected_files:
        if 'Vertex' in f:
            continue
        with reader(f) as ins:
            for l in ins:
                e = json.loads(l)  # edge
                if e['from'] == association_gid:
                    edges[f].append(e)
                if e['to'] == association_gid:
                    edges[f].append(e)
    normalized_edges = defaultdict(list)
    for f in normalized_files:
        with reader(f) as ins:
            for l in ins:
                e = json.loads(l)  # edge
                if e['from'] == association_gid:
                    normalized_edges[f].append(e)
                if e['to'] == association_gid:
                    normalized_edges[f].append(e)

    for k, v in edges.items():
        print(k, len(v))

    assert edges['outputs/pharmgkb/G2PAssociation_Genes_Gene.Edge.json.gz'][0]['to'] == 'ENSG00000188641', "Should be DPYD"
    assert edges['outputs/pharmgkb/G2PAssociation_Compounds_Compound.Edge.json.gz'][0]['to'] == 'Compound:TODO:fluorouracil', "Should be DPYD"
    assert sorted([e['to'] for e in edges['outputs/pharmgkb/G2PAssociation_Publications_Publication.Edge.json.gz']]) == ['Publication:ncbi.nlm.nih.gov/pubmed/10071185', 'Publication:ncbi.nlm.nih.gov/pubmed/24648345']

    assert len(edges['outputs/pharmgkb/G2PAssociation_GenomicFeatures_GenomicFeature.Edge.json.gz']) == 15, "Should have 15 edges to genomic feature"
    assert len([e for e in edges['outputs/pharmgkb/G2PAssociation_GenomicFeatures_GenomicFeature.Edge.json.gz'] if e['to'] == "GenomicFeature:b7fa874b108156b874c50502c5443b2ef9ee8b79"]) == 1, "Should have edge to DPYD*10"

    for k, v in normalized_edges.items():
        print(k, len(v))

    assert len(normalized_edges['outputs/compound/normalized.Compound_G2PAssociations_G2PAssociation.Edge.json.gz']) == len(edges['outputs/pharmgkb/Compound_G2PAssociations_G2PAssociation.Edge.json.gz']), 'Should have same number of normalized edges'
    assert len(normalized_edges['outputs/compound/normalized.G2PAssociation_Compounds_Compound.Edge.json.gz']) == len(edges['outputs/pharmgkb/G2PAssociation_Compounds_Compound.Edge.json.gz']), 'Should have same number of normalized edges'


def test_clinical_ann_metadata_981755803():
    """This association has a phenotype"""
    # pick specific association
    f = 'outputs/pharmgkb/G2PAssociation.Vertex.json.gz'
    association_gid = None
    with reader(f) as ins:
        for l in ins:
            a = json.loads(l)  # assoc
            cam = json.loads(a['data']['source_document'])  # clinical_ann_metadata
            # check it's content
            if cam['clinical_annotation_id'] == '981755803':
                association_gid = a['gid']
                assert a['data']['evidence_label'] == '1A'
                assert a['data']['source'] == 'pharmgkb'
                assert a['data']['source_url'] == 'https://www.pharmgkb.org/clinicalAnnotation/981755803'
    # check all edges
    edges = defaultdict(list)
    for f in expected_files:
        if 'Vertex' in f:
            continue
        with reader(f) as ins:
            for l in ins:
                e = json.loads(l)  # edge
                if e['from'] == association_gid:
                    edges[f].append(e)
                if e['to'] == association_gid:
                    edges[f].append(e)

    assert edges['outputs/pharmgkb/G2PAssociation_Genes_Gene.Edge.json.gz'][0]['to'] == 'ENSG00000001626', "Should be CFTR"
    assert edges['outputs/pharmgkb/G2PAssociation_Compounds_Compound.Edge.json.gz'][0]['to'] == 'Compound:TODO:ivacaftor', "Should be ivacaftor"
    assert sorted([e['to'] for e in edges['outputs/pharmgkb/G2PAssociation_Publications_Publication.Edge.json.gz']]) == \
        ['Publication:ncbi.nlm.nih.gov/pubmed/19846789',
         'Publication:ncbi.nlm.nih.gov/pubmed/21083385',
         'Publication:ncbi.nlm.nih.gov/pubmed/22047557',
         'Publication:ncbi.nlm.nih.gov/pubmed/22293084',
         'Publication:ncbi.nlm.nih.gov/pubmed/22942289',
         'Publication:ncbi.nlm.nih.gov/pubmed/23313410',
         'Publication:ncbi.nlm.nih.gov/pubmed/23590265',
         'Publication:ncbi.nlm.nih.gov/pubmed/23628510',
         'Publication:ncbi.nlm.nih.gov/pubmed/23757359',
         'Publication:ncbi.nlm.nih.gov/pubmed/23757361',
         'Publication:ncbi.nlm.nih.gov/pubmed/23891399',
         'Publication:ncbi.nlm.nih.gov/pubmed/24066763',
         'Publication:ncbi.nlm.nih.gov/pubmed/24461666',
         'Publication:ncbi.nlm.nih.gov/pubmed/25049054',
         'Publication:ncbi.nlm.nih.gov/pubmed/25145599',
         'Publication:ncbi.nlm.nih.gov/pubmed/25148434',
         'Publication:ncbi.nlm.nih.gov/pubmed/25171465',
         'Publication:ncbi.nlm.nih.gov/pubmed/25311995',
         'Publication:ncbi.nlm.nih.gov/pubmed/25473543',
         'Publication:ncbi.nlm.nih.gov/pubmed/25682022',
         'Publication:ncbi.nlm.nih.gov/pubmed/25755212',
         'Publication:ncbi.nlm.nih.gov/pubmed/26135562',
         'Publication:ncbi.nlm.nih.gov/pubmed/26568242',
         'Publication:ncbi.nlm.nih.gov/pubmed/27158673',
         'Publication:ncbi.nlm.nih.gov/pubmed/27745802',
         'Publication:ncbi.nlm.nih.gov/pubmed/27773592',
         'Publication:ncbi.nlm.nih.gov/pubmed/28611235',
         'Publication:ncbi.nlm.nih.gov/pubmed/28651844',
         'Publication:ncbi.nlm.nih.gov/pubmed/28711222']
    assert edges['outputs/pharmgkb/G2PAssociation_Phenotypes_Phenotype.Edge.json.gz'][0]['to'] == 'Phenotype:TODO:Cystic Fibrosis', "Should be Cystic Fibrosis"
    assert len(edges['outputs/pharmgkb/G2PAssociation_GenomicFeatures_GenomicFeature.Edge.json.gz']) == 1, "Should have 1 edge"
    assert len([e for e in edges['outputs/pharmgkb/G2PAssociation_GenomicFeatures_GenomicFeature.Edge.json.gz'] if e['to'] == "GenomicFeature:d924970d7daa4c04f790fee48e073dd30e4b954d"]) == 1, "Should have 1 edge to rs75527207"

    # check that all edges were normalized
    normalized_edges = defaultdict(list)
    for f in normalized_files:
        with reader(f) as ins:
            for l in ins:
                e = json.loads(l)  # edge
                if e['from'] == association_gid:
                    normalized_edges[f].append(e)
                if e['to'] == association_gid:
                    normalized_edges[f].append(e)

    for k, v in edges.items():
        print(k, len(v))

    for k, v in normalized_edges.items():
        print(k, len(v))

    assert len(normalized_edges['outputs/compound/normalized.Compound_G2PAssociations_G2PAssociation.Edge.json.gz']) == len(edges['outputs/pharmgkb/Compound_G2PAssociations_G2PAssociation.Edge.json.gz']), 'Should have same number of normalized edges'
    assert len(normalized_edges['outputs/compound/normalized.G2PAssociation_Compounds_Compound.Edge.json.gz']) == len(edges['outputs/pharmgkb/G2PAssociation_Compounds_Compound.Edge.json.gz']), 'Should have same number of normalized edges'
    assert len(normalized_edges['outputs/phenotype/normalized.Phenotype_G2PAssociations_G2PAssociation.Edge.json.gz']) == len(edges['outputs/pharmgkb/Phenotype_G2PAssociations_G2PAssociation.Edge.json.gz']), 'Should have same number of normalized edges'
    assert len(normalized_edges['outputs/phenotype/normalized.G2PAssociation_Phenotypes_Phenotype.Edge.json.gz']) == len(edges['outputs/pharmgkb/G2PAssociation_Phenotypes_Phenotype.Edge.json.gz']), 'Should have same number of normalized edges'


def _test_allele_edges():
    """As written, this test takes about 4 minutes to run.  Rename to run in development"""
    pharmgkb_allele_gids = []
    with reader('outputs/pharmgkb/G2PAssociation_Alleles_Allele.Edge.json.gz') as ins:
        for l in ins:
            e = json.loads(l)  # edge
            pharmgkb_allele_gids.append(e['to'])
    pharmgkb_allele_gids = list(set(pharmgkb_allele_gids))
    with reader('outputs/allele/Allele.Vertex.json.gz') as ins:
        for l in ins:
            v = json.loads(l)  # vertex
            if v['gid'] in pharmgkb_allele_gids:
                pharmgkb_allele_gids.remove(v['gid'])
            if len(pharmgkb_allele_gids) == 0:
                break
    print(len(pharmgkb_allele_gids))
    assert len(pharmgkb_allele_gids) == 0, "All pharmgkb alleles should be in normalized Alleles"
