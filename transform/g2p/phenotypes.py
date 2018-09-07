
from bmeg.vertex import Phenotype

# keep track of what we've already exported
EXPORTED_PHENOTYPES = []


def make_phenotype(term_id, term=None):
    """ return phenotype gid """
    return Phenotype(term_id=term_id, term=term)


def normalize(hit):
    """ return the hit modified replacing 'phenotypes'
    with phenotype_gids; phenotype_gids we haven't seen before """
    phenotypes = []
    association = hit['association']
    for phenotype in association.get('phenotypes', []):
        phenotypes.append(make_phenotype(phenotype['id'], phenotype.get('term', phenotype.get('description', None))))
    hit['phenotypes'] = phenotypes
    phenotype_gids = [p.gid() for p in phenotypes if p.gid() not in EXPORTED_PHENOTYPES]
    EXPORTED_PHENOTYPES.extend(phenotype_gids)
    return (hit, phenotype_gids)
