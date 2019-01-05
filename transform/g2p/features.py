
from bmeg.vertex import Allele, MinimalAllele
import bmeg.enrichers.gene_enricher as gene_enricher
from bmeg.vertex import Gene
import re

# keep track of what we've already exported
EXPORTED_ALLELES = []


def gene_gid(symbol):
    """ return gene gid """
    symbol = symbol.replace('Wild-Type', '').strip()
    gene = gene_enricher.get_gene(symbol)
    return Gene.make_gid(gene_id=gene['ensembl_gene_id'])


def allele(feature):
    """ return compound gid """
    params = {
        'genome': feature['referenceName'],
        'chromosome': feature['chromosome'],
        'start': feature['start'],
        'end': feature['end'],
        'reference_bases': feature.get('ref', None),
        'alternate_bases': feature.get('alt', None),
        'strand': '+',
        'hugo_symbol': feature.get('geneSymbol', None),
    }

    return Allele(**params)


def minimal_allele(feature):
    """ return compound gid """
    params = {
        'genome': feature.get('referenceName', None),
        'chromosome': feature.get('chromosome', None),
        'start': feature.get('start', None),
        'end': feature.get('end', None),
        'type': feature.get('biomarker_type', None),
        'name': feature.get('description', feature.get('name', None))
    }

    return MinimalAllele(**params)


def normalize(hit):
    """ return the hit modified replacing 'features'
    with allele_gids; allele_gids we haven't seen before """
    alleles = []
    missing_vertexes = []
    minimal_alleles = []
    allele_has_gene = []
    minimal_allele_has_gene = []
    for feature in hit['features']:
        if feature.get('provenance_rule', None) == 'gene_only':
            continue
        try:
            a = allele(feature)
            # skip if we already processed it
            if a.gid() in [_a.gid() for _a in alleles]:
                continue

            alleles.append(a)
            try:
                allele_has_gene.append((a.gid(), gene_gid(feature['geneSymbol'])))
            except Exception:
                missing_vertexes.append({'target_label': 'Gene', 'data': feature})
        except Exception as e:
            if 'start' in feature:
                print('>>could not create allele from {}'.format(str(e)))
            try:
                a = minimal_allele(feature)
                minimal_alleles.append(a)
                try:
                    # this check for geneSymbol should be in g2p, not here
                    description_parts = re.split(' +', feature['description'].strip())
                    geneSymbol = feature.get('geneSymbol', description_parts[0])
                    minimal_allele_has_gene.append((a.gid(), gene_gid(geneSymbol)))
                except Exception:
                    missing_vertexes.append({'target_label': 'Gene', 'data': feature})
            except Exception:
                missing_vertexes.append({'target_label': 'Allele', 'data': feature})

    hit['features'] = [a for a in alleles if a.gid() not in EXPORTED_ALLELES]
    hit['allele_has_gene'] = allele_has_gene
    hit['minimal_alleles'] = list(set([a for a in minimal_alleles if a.gid() not in EXPORTED_ALLELES]))
    hit['minimal_allele_has_gene'] = minimal_allele_has_gene

    allele_gids = [a.gid() for a in alleles]
    EXPORTED_ALLELES.extend(list(set(allele_gids)))
    minimal_allele_gids = [a.gid() for a in minimal_alleles]
    EXPORTED_ALLELES.extend(list(set(minimal_allele_gids)))
    return (hit, allele_gids, minimal_allele_gids, missing_vertexes)
