
import logging
from bmeg.vertex import Allele
from bmeg.util.logging import log_missing_vertex

# keep track of what we've already exported
EXPORTED_ALLELES = []


def allele_gid(feature):
    """ return compound gid """
    params = {
        'genome': feature['referenceName'],
        'chromosome': feature['chromosome'],
        'start': feature['start'],
        'end': feature['end'],
        'reference_bases': feature.get('ref', None),
        'alternate_bases': feature.get('alt', None),
    }

    return Allele.make_gid(**params)


def normalize(hit):
    """ return the hit modified replacing 'features'
    with allele_gids; allele_gids we haven't seen before """
    allele_gids = []
    for feature in hit['features']:
        try:
            allele_gids.append(allele_gid(feature))
        except Exception as e:
            logging.debug(e)
            log_missing_vertex({'label': 'Allele', 'feature': feature})

    hit['features'] = list(set(allele_gids))
    allele_gids = [gid for gid in allele_gids if gid not in EXPORTED_ALLELES]
    EXPORTED_ALLELES.extend(list(set(allele_gids)))
    return (hit, allele_gids)
