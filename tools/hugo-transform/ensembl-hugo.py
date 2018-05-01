#!/usr/bin/env python

# parses file found at ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc_complete_set.txt.gz
# and creates a map of ensembl to hugo symbols

# example usage -------------------

# python ensembl-hugo.py --hugo ~/Data/hugo/hgnc_complete_set.txt --out ~/Data/hugo/ensembl.json

import sys
import csv
import re
from pprint import pprint
import argparse
import json

def convert_line(state, line):
    if line['Status'] == 'Approved' and len(line['Ensembl Gene ID']) > 0:
        state[line['Ensembl Gene ID']] = line['Approved Symbol']

    return state

def convert_hugo(file):
    state = {}
    reader = csv.DictReader(file, delimiter='\t')
    for line in reader:
        state = convert_line(state, line)

    return state

def parse_args(args):
    args = args[1:]
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--hugo', type=str, help='path to the hugo source file')
    parser.add_argument('--out', type=str, help='Path to output file json')

    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    state = convert_hugo(open(options.hugo))
    json = json.dumps(state)
    with open(options.out, 'wb') as out:
        out.write(json)