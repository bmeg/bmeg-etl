#!/usr/bin/env python

import sys
import json
from bmeg.enrichers.drug_enricher import normalize

if __name__ == "__main__":
    out = normalize(sys.argv[1])
    if out is not None:
        out["_query"] = sys.argv[1]
        print(json.dumps(out))
