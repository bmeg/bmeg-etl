#!/usr/bin/python
# -*- coding: utf-8 -*-
import json
import os.path
import urllib.request

# load gene names
GENES = {}
ALIASES = {}

# get file if not already present
fname = 'non_alt_loci_set.json'
if not os.path.isfile(fname):
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/non_alt_loci_set.json'  # noqa
    urllib.request.urlretrieve(url, fname)

# trim payload, we only need symbol and ensembl
data = json.load(open(fname))
for doc in data['response']['docs']:
    gene = {
        'symbol': doc['symbol'],
        'ensembl_gene_id': doc.get('ensembl_gene_id', None),
        'entrez_id': doc.get('entrez_id', None)
        }
    GENES[doc['symbol']] = [gene]
    if gene['ensembl_gene_id']:
        if gene['ensembl_gene_id'] not in ALIASES:
            ALIASES[gene['ensembl_gene_id']] = []
        ALIASES[gene['ensembl_gene_id']].append(gene)
    if gene['entrez_id']:
        if gene['entrez_id'] not in ALIASES:
            ALIASES[gene['entrez_id']] = []
        ALIASES[gene['entrez_id']].append(gene)
    for alias in doc.get('alias_symbol', []):
        if alias not in ALIASES:
            ALIASES[alias] = []
        ALIASES[alias].append(gene)
    for prev in doc.get('prev_symbol', []):
        if prev not in ALIASES:
            ALIASES[prev] = []
        ALIASES[prev].append(gene)
data = None


def get_gene(identifier):
    """ return gene for identifier """
    genes = None
    for store in [GENES, ALIASES]:
        genes = store.get(identifier, None)
        if genes and len(genes) == 1:
            return genes
    if genes is None:
        raise ValueError(
            'gene reference does not exist {}'.format(identifier))  # noqa
    raise ValueError(
        'gene reference refers to multiple genes {}'.format(identifier))  # noqa
