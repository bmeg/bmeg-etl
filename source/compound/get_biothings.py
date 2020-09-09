#!/usr/bin/env python

import ujson as json
import os
import time
from glob import glob
from bmeg.enrichers.drug_enricher import normalize_biothings

def retry(count, f, *args, **kwargs ):
    while count > 0:
        out = f(*args, **kwargs)
        if out is not None:
            return out
        count -= 1
        time.sleep(1)
    raise Exception("Retry Failed")

def run_biothings(
                    table_dir="./reference/compound/",
                    outfile="./source/compound/biothings.json"):
    ids = {}
    with open(outfile, "w") as ohandle:
        for f in glob(os.path.join(table_dir, "*.table")):
            with open(f) as handle:
                for line in handle:
                    row = line.rstrip().split("\t")
                    id = row[1]
                    if id not in ids:
                        out = normalize_biothings(id)
                        if out is None:
                            print("Error finding BioThing Entry: %s" %id)
                        else:
                            ohandle.write("%s\n" % (json.dumps(out)))
                            ids[id] = True

if __name__ == "__main__":
    run_biothings()
