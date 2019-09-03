import pydash

from bmeg import Phenotype, Project
from bmeg.enrichers.phenotype_enricher import phenotype_factory


# keep track of what we've already exported
EXPORTED_PHENOTYPES = {}


def make_phenotype(term_id, term=None):
    """ return phenotype gid """
    return Phenotype(term_id=term_id,
                     term=term,
                     id=Phenotype.make_gid(term_id),
                     project_id=Project.make_gid('Reference'))


def normalize(hit):
    """ return the hit modified replacing 'phenotypes'
    with phenotype_gids; phenotype_gids we haven't seen before """
    phenotypes = set([])
    # we just want the main phenotype not the parent terms
    if 'molecularmatch' in hit['source']:
        for tag in pydash.get(hit, 'molecularmatch_trials.tags', pydash.get(hit, 'molecularmatch.tags', [])):
            if tag.get('filterType', '') == 'include' and tag.get('priority', 0) == 1 and tag.get('facet', '') == 'CONDITION':
                phenotypes.add(phenotype_factory(tag['term']))
    # handle special case for pmkb cancer genes
    elif all([hit['source'] == 'pmkb',
              len(pydash.get(hit, 'association.phenotypes', [])) == 148,
              pydash.get(hit, 'association.description', '') == 'This gene is a known cancer gene.']):
        phenotypes.add(phenotype_factory('cancer'))
        hit['association']['oncogenic'] = 'known cancer gene'
    # general case
    else:
        association = hit['association']
        for phenotype in association.get('phenotypes', []):
            if phenotype['id'].startswith('MONDO'):
                phenotypes.add(make_phenotype(phenotype['id'], phenotype.get('term', phenotype.get('description', None))))
            else:
                # should we use id or term here?
                phenotypes.add(phenotype_factory(phenotype['id']))

    hit['phenotypes'] = []
    for p in phenotypes:
        if p.gid() in EXPORTED_PHENOTYPES:
            continue
        hit['phenotypes'].append(p)
        EXPORTED_PHENOTYPES[p.gid()] = True
    phenotype_gids = list(set([p.gid() for p in phenotypes]))
    return (hit, phenotype_gids)
