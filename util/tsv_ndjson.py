#!/usr/bin/env python

import sys
import csv
import json

with open(sys.argv[1], "rt") as handle:
    header = None
    reader = csv.reader(handle, delimiter='\t', quotechar='"')
    for row in reader:
        if header is None:
            header = row
        else:
            out = dict(zip(header, row))
            print(json.dumps(out))