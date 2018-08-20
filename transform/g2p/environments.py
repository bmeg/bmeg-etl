
from bmeg.vertex import Compound

# keep track of what we've already exported
EXPORTED_COMPOUNDS = []


def compound(term_id, term=None):
    """ return compound gid """
    return Compound(term_id=term_id)


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with compound_gids; compound_gids we haven't seen before """
    compounds = []
    for environment in hit.get('environments', []):
        compounds.append(compound(environment))
    hit['environments'] = compounds
    compound_gids = [compound.gid() for compound in compounds if compound.gid() not in EXPORTED_COMPOUNDS]
    EXPORTED_COMPOUNDS.extend(compound_gids)
    return (hit, compound_gids)
