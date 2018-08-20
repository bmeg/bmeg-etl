
from bmeg.vertex import Publication

# keep track of what we've already exported
EXPORTED_PUBLICATIONS = []


def publication(url, description=None):
    """ return publication gid """
    return Publication(url=url, description=description)


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with publication_gids; publication_gids we haven't seen before """
    publications = []
    # format evidence as bmeg friendly
    association = hit['association']
    evidence = association['evidence'][0]
    if evidence.get('info', None):
        for url in evidence['info'].get('publications', []):
            publications.append(publication(url))
    hit['publications'] = publications
    publication_gids = [p.gid() for p in publications if p.gid() not in EXPORTED_PUBLICATIONS]
    EXPORTED_PUBLICATIONS.extend(publication_gids)
    return (hit, publication_gids)
