
from bmeg.vertex import Compound
from bmeg.enrichers.drug_enricher import compound_factory
# keep track of what we've already exported
EXPORTED_COMPOUNDS = []


def compound(environment):
    """ return compound gid """
    if 'term' in environment and environment.get('id', None):
        return Compound(name=environment['description'], term=environment['term'], term_id=environment['id'])
    else:
        return compound_factory(name=environment['description'])


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with compound_gids; compound_gids we haven't seen before """
    compounds = set([])
    association = hit['association']
    for environment in association.get('environmentalContexts', []):
        compounds.add(compound(environment))
    hit['environments'] = [compound for compound in compounds if compound.gid() not in EXPORTED_COMPOUNDS]
    compound_gids = [compound.gid() for compound in compounds]
    EXPORTED_COMPOUNDS.extend(compound_gids)
    return (hit, compound_gids)
