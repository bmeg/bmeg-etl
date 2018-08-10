#!/usr/bin/env python3

import argparse
import pandas
import csv
import gzip

from bmeg.vertex import GeneExpression
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--matrix', type=str, help='')
    parser.add_argument("--output")

    args = parser.parse_args()

    emitter = JSONEmitter(args.output)

    #df = pandas.read_csv(args.matrix, sep="\t", index_col=0, dtype="a", engine="c", iterator=True, compression="gzip")
    #df = pandas.read_csv(args.matrix, sep="\t", index_col=0, iterator=True, compression="gzip")
    print("Done Loading")
    with gzip.open(args.matrix, "rt") as handle:
        df = csv.reader(handle, delimiter="\t")
        samples = next(df)[1:]
        for k in df:
            assert(len(samples) == len(values))
            #print(gene, dict(zip(samples, values)))
