
from bmeg import Compound, Project
from bmeg.enrichers.drug_enricher import compound_factory

# keep track of what we've already exported
EXPORTED_COMPOUNDS = []


def compound(environment):
    """ return compound gid """
    if 'term' in environment and environment.get('id', None):
        return Compound(
            name=environment['description'], term=environment['term'], term_id=environment['id'],
            id=Compound.make_gid(environment['id']),
            project_id=Project.make_gid("Reference"),
        )
    else:
        return compound_factory(name=environment['description'])


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with compound_gids; compound_gids we haven't seen before """
    compounds = set([])
    association = hit['association']
    for environment in association.get('environmentalContexts', []):
        compounds.add(compound(environment))

    hit['environments'] = []
    dups = []
    for c in compounds:
        if c.gid() in EXPORTED_COMPOUNDS:
            continue
        if c.gid() in dups:
            continue
        hit['environments'].append(c)
        dups.append(c.gid())
    compound_gids = [compound.gid() for compound in compounds]
    EXPORTED_COMPOUNDS.extend(compound_gids)
    return (hit, compound_gids)
