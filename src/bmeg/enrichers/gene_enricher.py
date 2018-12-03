""" lookup gene alias and id """
import json
import os.path
import urllib.request

# load gene names
GENES = {}
ALIASES = {}


def append(key, val, collection):
    if not key:
        return
    if key not in collection:
        collection[key] = [val]
        return
    found = False
    for item in collection[key]:
        if item == val:
            found = True
    if not found:
        collection[key].append(val)
    return


# get file if not already present
fname = 'source/gene_enricher/hgnc_complete_set.json'
if not os.path.isfile(fname):
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json'
    urllib.request.urlretrieve(url, fname)

# trim payload, we only need symbol and ensembl
data = json.load(open(fname))
for doc in data['response']['docs']:
    gene = {
        'symbol': doc['symbol'],
        'ensembl_gene_id': doc.get('ensembl_gene_id', None),
        'entrez_id': doc.get('entrez_id', None),
        'hgnc_id': doc.get('hgnc_id', None)
    }

    append(doc['symbol'], gene, GENES)
    if gene['ensembl_gene_id']:
        append(doc['ensembl_gene_id'], gene, ALIASES)
    if gene['entrez_id']:
        append(doc['entrez_id'], gene, ALIASES)
    if gene['hgnc_id']:
        append(doc['hgnc_id'], gene, ALIASES)
    for alias in doc.get('alias_symbol', []):
        append(alias, gene, ALIASES)
    for prev in doc.get('prev_symbol', []):
        append(prev, gene, ALIASES)
data = None


def get_gene(identifier):
    """ return gene for identifier """
    genes = None
    for store in [GENES, ALIASES]:
        genes = store.get(identifier, None)
        if genes and len(genes) == 1:
            return genes[0]
    if genes is None:
        raise ValueError(
            "gene reference does not exist: '{}'".format(identifier))  # noqa
    raise ValueError(
        "gene reference refers to multiple genes: '{}'".format(identifier))  # noqa
