#!/usr/bin/env python3

import os
import logging

import bmeg.requests
from bmeg.util.cli import default_argument_parser


client = bmeg.requests.Client("pfam")
parser = default_argument_parser()
args = parser.parse_args()

output_dir = "source/pfam/xmls"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists("source/pfam/id_list.txt"):
    raise Exception("missing ID list. Run transform/pfam/list.py")

ids = open("source/pfam/id_list.txt").read().splitlines()

for i in ids:
    url = "http://pfam.xfam.org/family?output=xml&acc=%s" % i
    logging.debug("downloading {}", url)
    handle = client.get(url)
    txt = handle.text
    f = os.path.join(output_dir, "%s.xml" % (i))
    with open(f, "w") as handle:
        handle.write(txt)
