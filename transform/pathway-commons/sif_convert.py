#!/usr/bin/env python

from __future__ import print_function

import argparse
import csv
import gzip
import json
import logging
import sys

from bmeg.pathway_pb2 import ProteinInteraction
from google.protobuf import json_format


def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    return json.dumps(msg)


def parse_gene_map(path):
    name_table = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if len(row['Ensembl Gene ID']) > 0:
                name_table[row['Approved Symbol']] = row['Ensembl Gene ID']
                if len(row['Previous Symbols']) > 0:
                    for ps in row['Previous Symbols'].split("|"):
                        name_table[ps] = row['Ensembl Gene ID']
    return name_table


def parse_sif(sif, name_table):
    with gzip.GzipFile(sif) as handle:
        for line in handle:
            row = line.rstrip().split("\t")
            if not (row[0].startswith("CHEBI:") or row[2].startswith("CHEBI:")
                    or row[2].startswith("PubChem:")):
                pi = ProteinInteraction()
                if row[0] not in name_table:
                    logging.warning("%s not found" % (row[0]))
                elif row[2] not in name_table:
                    logging.warning("%s not found" % (row[2]))
                else:
                    pi.src = name_table[row[0]]
                    pi.dst = name_table[row[2]]
                    pi.interaction = row[1]
                    print(message_to_json(pi))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser()
    # http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons9.All.hgnc.sif.gz
    parser.add_argument("-s", "--sif", required=True)
    # ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
    parser.add_argument("-g", "--gene-map", required=True)
    parser.add_argument("-o", "--output", default=sys.stdout)
    args = parser.parse_args()

    if args.output != sys.stdout:
        sys.stdout = open(args.output, 'w')

    name_table = parse_gene_map(args.gene_map)
    parse_sif(args.sif, name_table)
