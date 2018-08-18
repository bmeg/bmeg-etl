
from bmeg.vertex import Publication

# keep track of what we've already exported
EXPORTED_PUBLICATIONS = []


def publication_gid(url):
    """ return publication gid """
    return Publication.make_gid(url=url)


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with publication_gids; publication_gids we haven't seen before """
    publication_gids = []
    # format evidence as bmeg friendly
    association = hit['association']
    evidence = association['evidence'][0]
    if evidence.get('info', None):
        for url in evidence['info'].get('publications', []):
            publication_gids.append(publication_gid(url))
    hit['publications'] = publication_gids
    publication_gids = [gid for gid in publication_gids if gid not in EXPORTED_PUBLICATIONS]
    EXPORTED_PUBLICATIONS.extend(publication_gids)
    return (hit, publication_gids)
