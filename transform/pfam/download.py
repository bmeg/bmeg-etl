#!/usr/bin/env python3

import os
import sys
import logging
import requests
try:
    from tqdm import tqdm
except ImportError:
    tdqm = lambda x:x

logging.basicConfig(level=logging.INFO)

output_dir = sys.argv[1]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

ids = []
handle = requests.get("http://pfam.xfam.org/families?output=text")
for line in handle.iter_lines():
    line = line.decode()
    if not line.startswith("#"):
        row = line.split("\t")
        if len(row[0]) > 1:
            ids.append(row[0])

for i in tqdm(ids):
    f = os.path.join(output_dir, "%s.xml" % (i))
    if not os.path.exists(f):
        url = "http://pfam.xfam.org/family?output=xml&acc=%s" % i
        logging.debug("downloading %s to %s" % (url, f))
        handle = requests.get(url)
        txt = handle.text
        with open(f, "w") as handle:
            handle.write(txt)
    else:
        logging.debug("skipping %s" % (i))
