import logging
import re

from transform.g2p.genes import gene_gid
from bmeg import Allele, GenomicFeature, Project

# keep track of what we've already exported
EXPORTED_ALLELES = {}
ALLELE_HAS_GENE_CACHE = {}
GENOMIC_FEATURE_HAS_GENE_CACHE = {}


def allele(feature):
    """ return alle """
    params = {
        'genome': feature['referenceName'],
        'chromosome': feature['chromosome'],
        'start': feature['start'],
        'reference_bases': feature['ref'],
        'alternate_bases': feature['alt'],
        'id': Allele.make_gid(
            feature['referenceName'], feature['chromosome'],
            feature['start'],
            feature['ref'], feature['alt']
        ),
        'project_id': Project.make_gid("Reference")
    }
    a = Allele(**params)
    a.validate()
    return a


def genomic_feature(feature):
    """ return genomic feature """
    params = {
        'genome': feature.get('referenceName', None),
        'chromosome': feature.get('chromosome', None),
        'start': feature.get('start', None),
        'end': feature.get('end', None),
        'type': feature.get('biomarker_type', None),
        'name': feature.get('description', feature.get('name', None)),
        'id': GenomicFeature.make_gid(
            feature.get('referenceName', None), feature.get('chromosome', None),
            feature.get('start', None), feature.get('end', None),
            feature.get('biomarker_type', None), feature.get('name', None)
        ),
        'project_id': Project.make_gid("Reference")
    }
    gf = GenomicFeature(**params)
    gf.validate()
    return gf


def normalize(hit):
    """ return the hit modified replacing 'features'
    with allele_gids; allele_gids we haven't seen before """
    alleles = set([])
    missing_vertexes = []
    genomic_features = set([])
    allele_has_gene = set([])
    genomic_feature_has_gene = set([])
    for feature in hit['features']:
        if feature.get('provenance_rule', None) == 'gene_only':
            continue
        try:
            a = allele(feature)
            alleles.add(a)
        except Exception:
            logging.debug("unable to convert feature into an Allele; creating a GenomicFeature")
            try:
                a = genomic_feature(feature)
                genomic_features.add(a)
                try:
                    # this check for geneSymbol should be in g2p, not here
                    description_parts = re.split(' +', feature['description'].strip())
                    geneSymbol = feature.get('geneSymbol', description_parts[0])
                    genomic_feature_has_gene.add((a.gid(), gene_gid(geneSymbol)))
                except Exception:
                    missing_vertexes.append({'target_label': 'Gene', 'data': feature})
            except Exception:
                missing_vertexes.append({'target_label': 'Allele', 'data': feature})

    hit['features'] = []
    for a in alleles:
        if a.gid() not in EXPORTED_ALLELES:
            hit['features'].append(a)
            EXPORTED_ALLELES[a.gid()] = True

    hit['genomic_features'] = []
    for a in genomic_features:
        if a.gid() not in EXPORTED_ALLELES:
            hit['genomic_features'].append(a)
            EXPORTED_ALLELES[a.gid()] = True

    hit['allele_has_gene'] = []
    for a in allele_has_gene:
        if a not in ALLELE_HAS_GENE_CACHE:
            hit['allele_has_gene'].append(a)
            ALLELE_HAS_GENE_CACHE[a] = True

    hit['genomic_feature_has_gene'] = []
    for a in genomic_feature_has_gene:
        if a not in GENOMIC_FEATURE_HAS_GENE_CACHE:
            hit['genomic_feature_has_gene'].append(a)
            GENOMIC_FEATURE_HAS_GENE_CACHE[a] = True

    allele_gids = list(set([a.gid() for a in alleles]))
    genomic_feature_gids = list(set([a.gid() for a in genomic_features]))

    return (hit, allele_gids, genomic_feature_gids, missing_vertexes)
