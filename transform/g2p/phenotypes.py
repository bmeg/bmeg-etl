
from bmeg import Phenotype, Project

# keep track of what we've already exported
EXPORTED_PHENOTYPES = []


def make_phenotype(term_id, term=None):
    """ return phenotype gid """
    return Phenotype(term_id=term_id, term=term,
                     id=Phenotype.make_gid(term_id),
                     project_id=Project.make_gid("Reference"))


def normalize(hit):
    """ return the hit modified replacing 'phenotypes'
    with phenotype_gids; phenotype_gids we haven't seen before """
    phenotypes = set([])
    association = hit['association']
    for phenotype in association.get('phenotypes', []):
        phenotypes.add(make_phenotype(phenotype['id'], phenotype.get('term', phenotype.get('description', None))))
    hit['phenotypes'] = []
    dups = []
    for p in phenotypes:
        if p.gid() in EXPORTED_PHENOTYPES:
            continue
        if p.gid() in dups:
            continue
        hit['phenotypes'].append(p)
        dups.append(p.gid())
    phenotype_gids = list(set([p.gid() for p in phenotypes]))
    EXPORTED_PHENOTYPES.extend(phenotype_gids)
    return (hit, phenotype_gids)
