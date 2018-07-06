#!/usr/bin/env python

import argparse
import gzip
import json
import re
import sys

from google.protobuf import json_format
import bmeg.rna_pb2 as rna

parser = argparse.ArgumentParser(
    prog="convert-expression",
    description="Convert a tab-separated expression matrix into a stream of GeneExpression JSON messages.")

parser.add_argument("source",
    help="Describes the source of the data e.g. 'ccle'.")

parser.add_argument("scale",
    choices=["unknown", "read_count", "tpkm", "rpkm", "fpkm"],
    help="Describes the scale of the expression data e.g. 'rpkm'.")

parser.add_argument("input",
    help="The input TSV file to convert.")

parser.add_argument("output",
    help="The output file to write to.")

args = parser.parse_args()

opener = open
if args.input.endswith(".gz"):
    opener = gzip.open
    
# Convert string arg to rna enum value
scale = {
    "unknown": rna.UNKNOWN,
    "read_count": rna.READ_COUNT,
    "tpkm": rna.TPKM,
    "rpkm": rna.RPKM,
    "fpkm": rna.FPKM,
}[args.scale]

print "loading table"

rows = []
fh = opener(args.input)
fh.readline() # skip version number

# Table dimensions header
dim = fh.readline().strip().split("\t")
num_rows, num_cols = int(dim[0]), int(dim[1])

header = fh.readline() # header with sample names
header = header.strip().split("\t")
samples = header[2:]

# Begin expression matrix loading.
#
# Rows are genes, columns are samples.
# Build one sample message at a time.
# For each sample column, loop over all rows, gather the expression values,
# and dump the message.
#
# The data is large, so this avoids loading it all into memory
# by looping over the file multiple times using seek/tell.
start = fh.tell()

with open(args.output, "w") as out:
    for sample_i, sample in enumerate(samples):
        print "processing sample", sample
        fh.seek(start)

        msg = rna.GeneExpression()
        msg.id = sample
        # msg.id = args.source + "_expression:" + sample
        msg.source = args.source
        msg.biosample_id = sample
        msg.scale = scale

        for line in fh:
            line = line.strip().split("\t")
            name, desc, exprs = line[0], line[1], line[2:]
            if name != "":
                msg.expressions[name] = float(exprs[sample_i])

        s = json.dumps(json_format.MessageToDict(msg))
        out.write(s + "\n")
