#!/usr/bin/env python

import os
import sys
import json
import gzip

o = list()

for root, dirs, files in os.walk(sys.argv[1]):
    for name in files:
        if "Compound" in name:
            f = os.path.join(root, name)
            with gzip.open(f) as handle:
                for line in handle:
                    d = json.loads(line)
                    if "to" in d:
                        x = d.get("to", "")
                        if x.startswith("Compound"):
                            o.append(x)
                        x = d.get("from", "")
                        if x.startswith("Compound"):
                            o.append(x)
                    else:
                        x = d.get("gid", "")
                        if x.startswith("Compound"):
                            o.append(x)

for i in set(o):
    print(i.replace("Compound:", ""))
