#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
import struct
from array import array

from bmeg.vertex import GeneExpression
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter


def matrix_parse(handle, emitter, args):
    df = csv.reader(handle, delimiter="\t")
    samples = next(df)[1:]
    genes = []
    #write values out to temp matrix
    t = open("tmp_matrix.data", "wb")
    for k in df:
        gene = k[0]
        values = k[1:]
        assert(len(samples) == len(values))
        v = array("d", (float(a) for a in values))
        v.tofile(t)
        genes.append(gene)
    t.close()
    #read back, steping through the rows to do a transpose
    t = open("tmp_matrix.data", "rb")
    rowSize = 8 * len(samples)
    for i, s in enumerate(samples):
        values = []
        for j, g in enumerate(genes):
            t.seek(i*8 + j*rowSize)
            o = struct.unpack('d', t.read(8))
            values.append(o[0])
        out = dict(zip(genes, values))
        emitter.emit_vertex(GeneExpression(
            id=s,
            source=args.source,
            scale=args.scale,
            method=args.method,
            values=out
        ))
    t.close()
    os.unlink("tmp_matrix.data")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--matrix', type=str, help='')
    parser.add_argument("--output", default="rna")
    parser.add_argument("--scale", default="TPM")
    parser.add_argument("--method", default="Illumina HiSeq")
    parser.add_argument("--source", default="NA")
    parser.add_argument("--gz", action="store_true", default=False)

    args = parser.parse_args()

    emitter = JSONEmitter(args.output)
    if args.gz:
        with gzip.open(args.matrix, "rt") as handle:
            matrix_parse(handle, emitter, args)
    else:
        with open(args.matrix, "rt") as handle:
            matrix_parse(handle, emitter, args)
