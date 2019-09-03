
from bmeg import Publication, Project


EXPORTED_PUBLICATIONS = {}


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with publication_gids; publication_gids we haven't seen before """
    publications = set([])
    # format evidence as bmeg friendly
    association = hit['association']
    evidence = association['evidence'][0]
    if evidence.get('info', None):
        for url in evidence['info'].get('publications', []):
            publications.add(
                Publication(url=url.strip(), title=None, abstract=None, text=None, date=None, author=None, citation=None,
                            id=Publication.make_gid(url.strip()), project_id=Project.make_gid("Reference"))
            )

    # for p in publications:
    #     if p.gid() not in EXPORTED_PUBLICATIONS:
    #         hit['publications'].append(p)
    #         EXPORTED_PUBLICATIONS[p.gid()] = True

    publication_gids = list(set([p.gid() for p in publications]))
    return (hit, publication_gids)
