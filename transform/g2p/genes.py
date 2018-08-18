
import bmeg.maf.gene_enricher as gene_enricher
from bmeg.vertex import Gene
from bmeg.util.logging import log_missing_vertex
import logging
import json

EXPORTED_GENES = []


def gene_gid(symbol):
    """ return gene gid """
    symbol = symbol.replace('Wild-Type', '').strip()
    genes = gene_enricher.get_gene(symbol)
    return Gene.make_gid(gene_id=genes[0]['ensembl_gene_id'])


def normalize(hit):
    """ return the hit modified replacing 'genes' with create_genes;
    gene gids we haven't seen before """
    gene_gids = []
    for symbol in hit['genes']:
        try:
            gene_gids.append(gene_gid(symbol))
        except Exception as e:
            log_missing_vertex({'label': 'Gene', 'symbol': symbol})
    hit['genes'] = gene_gids
    gene_gids = [gid for gid in gene_gids if gid not in EXPORTED_GENES]
    EXPORTED_GENES.extend(gene_gids)
    return (hit, gene_gids)
