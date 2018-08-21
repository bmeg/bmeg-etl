
import logging
from bmeg.vertex import Allele

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


def normalize(hit):
    """ return the hit modified replacing 'features'
    with allele_gids; allele_gids we haven't seen before """
    alleles = []
    missing_vertexes = []
    for feature in hit['features']:
        if feature.get('provenance_rule', None) == 'gene_only':
            continue
        try:
            alleles.append(allele(feature))
        except Exception as e:
            logging.debug(e)
            missing_vertexes.append({'target_label': 'Allele', 'data': feature})

    hit['features'] = alleles
    allele_gids = [a.gid() for a in alleles if a.gid() not in EXPORTED_ALLELES]
    EXPORTED_ALLELES.extend(list(set(allele_gids)))
    return (hit, allele_gids, missing_vertexes)
