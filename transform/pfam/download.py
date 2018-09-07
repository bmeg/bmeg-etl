#!/usr/bin/env python3

import os
import tarfile

import bmeg.requests
from bmeg.util.cli import default_argument_parser


client = bmeg.requests.Client("pfam")
parser = default_argument_parser()
parser.add_argument(
    "--archive",
    "-A",
    action="store_true",
    default=False,
    help="create a gzipped TAR archive of all the output files")

args = parser.parse_args()

output_dir = "outputs/pfam"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists("outputs/pfam/id_list.txt"):
    raise Exception("missing ID list. Run transform/pfam/list.py")

ids = open("outputs/pfam/id_list.txt").read().splitlines()

if args.archive:
    tar = tarfile.open(os.path.join(output_dir, "pfam.tar.gz"), "w:gz")

for i in ids:
    handle = client.get("http://pfam.xfam.org/family?output=xml&acc=%s" % i)
    txt = handle.text
    f = os.path.join(output_dir, "%s.xml" % (i))
    with open(f, "w") as handle:
        handle.write(txt)

    if args.archive:
        tar.add(f, arcname=os.path.basename(f))

if args.archive:
    tar.close()
