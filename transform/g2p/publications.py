
from bmeg import Publication, Project

# keep track of what we've already exported
EXPORTED_PUBLICATIONS = []


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with publication_gids; publication_gids we haven't seen before """
    publications = []
    # format evidence as bmeg friendly
    association = hit['association']
    evidence = association['evidence'][0]
    if evidence.get('info', None):
        for url in evidence['info'].get('publications', []):
            publications.append(
                Publication(url=url.strip(), title=None, abstract=None, text=None, date=None, author=None, citation=None,
                            id=Publication.make_gid(url.strip()), project_id=Project.make_gid("Reference"))
            )
    hit['publications'] = publications
    publication_gids = [p.gid() for p in publications if p.gid() not in EXPORTED_PUBLICATIONS]
    EXPORTED_PUBLICATIONS.extend(publication_gids)
    return (hit, publication_gids)
