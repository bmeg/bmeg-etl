
import bmeg.enrichers.gene_enricher as gene_enricher
from bmeg.vertex import Gene
import logging

EXPORTED_GENES = []


def gene_gid(symbol):
    """ return gene gid """
    symbol = symbol.replace('Wild-Type', '').strip()
    genes = gene_enricher.get_gene(symbol)
    return Gene.make_gid(gene_id=genes[0]['ensembl_gene_id'])


def normalize(hit):
    """ return the hit modified replacing 'genes' with create_genes;
    gene gids we haven't seen before """
    gene_gids = set([])
    missing_vertexes = []
    for symbol in hit['genes']:
        try:
            gene_gids.add(gene_gid(symbol))
        except Exception as e:
            logging.debug(e)
            missing_vertexes.append({'target_label': 'Gene', 'data': {'symbol': symbol}})
    for feature in hit.get('features', []):
        if feature.get('provenance_rule', None) != 'gene_only':
            continue
        try:
            gene_gids.add(gene_gid(feature['geneSymbol']))
        except Exception as e:
            logging.debug(e)
            missing_vertexes.append({'target_label': 'Gene', 'data': feature})

    hit['genes'] = gene_gids
    gene_gids = [gid for gid in gene_gids if gid not in EXPORTED_GENES]
    EXPORTED_GENES.extend(gene_gids)
    return (hit, gene_gids, missing_vertexes)
