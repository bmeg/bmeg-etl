from bmeg.enrichers.drug_enricher import compound_factory

# keep track of what we've already exported
EXPORTED_COMPOUNDS = {}


def make_compound(environment):
    """ return compound gid """
    if environment.get('id', None):
        return compound_factory(name=environment['id'])
    else:
        return compound_factory(name=environment['description'])


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with compound_gids; compound_gids we haven't seen before """
    compounds = set([])
    association = hit['association']
    for environment in association.get('environmentalContexts', []):
        compounds.add(make_compound(environment))

    hit['environments'] = []
    for c in compounds:
        if c.gid() in EXPORTED_COMPOUNDS:
            continue
        hit['environments'].append(c)
        EXPORTED_COMPOUNDS[c.gid()] = True
    compound_gids = list(set([compound.gid() for compound in compounds]))
    return (hit, compound_gids)
