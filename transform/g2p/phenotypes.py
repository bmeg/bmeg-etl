
from bmeg.vertex import Phenotype

# keep track of what we've already exported
EXPORTED_PHENOTYPES = []


def make_phenotype(term_id, term=None):
    """ return phenotype gid """
    return Phenotype(term_id=term_id, term=term)


def normalize(hit):
    """ return the hit modified replacing 'phenotypes'
    with phenotype_gids; phenotype_gids we haven't seen before """
    phenotypes = set([])
    association = hit['association']
    for phenotype in association.get('phenotypes', []):
        phenotypes.add(make_phenotype(phenotype['id'], phenotype.get('term', phenotype.get('description', None))))
    hit['phenotypes'] = [p for p in phenotypes if p.gid() not in EXPORTED_PHENOTYPES]
    phenotype_gids = [p.gid() for p in phenotypes]
    EXPORTED_PHENOTYPES.extend(phenotype_gids)
    return (hit, phenotype_gids)
