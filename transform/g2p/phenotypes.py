
from bmeg.vertex import Phenotype

# keep track of what we've already exported
EXPORTED_PHENOTYPES = []


def phenotype_gid(term_id, term=None):
    """ return phenotype gid """
    return Phenotype.make_gid(term_id=term_id)


def normalize(hit):
    """ return the hit modified replacing 'phenotypes'
    with phenotype_gids; phenotype_gids we haven't seen before """
    phenotype_gids = []
    association = hit['association']
    for phenotype in association.get('phenotypes', []):
        phenotype_gids.append(phenotype_gid(phenotype['id']))
    hit['phenotypes'] = phenotype_gids
    phenotype_gids = [gid for gid in phenotype_gids if gid not in EXPORTED_PHENOTYPES]
    EXPORTED_PHENOTYPES.extend(phenotype_gids)
    return (hit, phenotype_gids)
