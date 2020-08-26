#!/usr/bin/env python

import sys
import ujson as json
from bmeg.enrichers.drug_enricher import normalize, search_pubchem, normalize_biothings


def process(name):
    cid = search_pubchem(name)
    compound = None
    if cid:
        compound = normalize_biothings(cid)
    if compound is None:
        compound = normalize_biothings(name, fuzzy=False)
    if compound is not None:
        print("%s\t%s\t%s" % (name, compound["id"], json.dumps(compound)))
    else:
        print("%s\t\t" % (name))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        for n in sys.argv[1:]:
            process(n)
    else:
        for line in sys.stdin:
            n = line.rstrip()
            process(n)
