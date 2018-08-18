
from bmeg.vertex import Compound

# keep track of what we've already exported
EXPORTED_COMPOUNDS = []


def compound_gid(term_id, term=None):
    """ return compound gid """
    return Compound.make_gid(term_id=term_id)


def normalize(hit):
    """ return the hit modified replacing 'environmentalContexts'
    with compound_gids; compound_gids we haven't seen before """
    compound_gids = []
    for environmentalContext in hit.get('environmentalContexts', []):
        compound_gids.append(compound_gid(environmentalContext))
    hit['environmentalContexts'] = compound_gids
    compound_gids = [gid for gid in compound_gids if gid not in EXPORTED_COMPOUNDS]
    EXPORTED_COMPOUNDS.extend(compound_gids)
    return (hit, compound_gids)
