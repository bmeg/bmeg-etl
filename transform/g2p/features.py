
import logging
from bmeg.vertex import Allele, MinimalAllele

# keep track of what we've already exported
EXPORTED_ALLELES = []


def allele(feature):
    """ return compound gid """
    params = {
        'genome': feature['referenceName'],
        'chromosome': feature['chromosome'],
        'start': feature['start'],
        'end': feature['end'],
        'reference_bases': feature.get('ref', None),
        'alternate_bases': feature.get('alt', None),
    }

    return Allele(**params)


def minimal_allele(feature):
    """ return compound gid """
    params = {
        'genome': feature['referenceName'],
        'chromosome': feature['chromosome'],
        'start': feature.get('start', None),
        'end': feature.get('end', None),
        'type': feature.get('biomarker_type', None),
        'name': feature.get('description', None)
    }

    return MinimalAllele(**params)


def normalize(hit):
    """ return the hit modified replacing 'features'
    with allele_gids; allele_gids we haven't seen before """
    alleles = []
    missing_vertexes = []
    minimal_alleles = []
    for feature in hit['features']:
        if feature.get('provenance_rule', None) == 'gene_only':
            continue
        try:
            alleles.append(allele(feature))
        except Exception:
            try:
                minimal_alleles.append(minimal_allele(feature))
            except Exception as e:
                logging.exception(e)
                missing_vertexes.append({'target_label': 'Allele', 'data': feature})

    hit['features'] = alleles
    hit['minimal_alleles'] = minimal_alleles
    allele_gids = [a.gid() for a in alleles if a.gid() not in EXPORTED_ALLELES]
    EXPORTED_ALLELES.extend(list(set(allele_gids)))
    minimal_allele_gids = [a.gid() for a in minimal_alleles if a.gid() not in EXPORTED_ALLELES]
    EXPORTED_ALLELES.extend(list(set(minimal_allele_gids)))

    return (hit, allele_gids, minimal_allele_gids, missing_vertexes)
