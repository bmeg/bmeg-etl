#!/usr/bin/env python

import sys
import json

def parse_table(path):
    table = {}
    with open(path) as compounds:
        for line in compounds:
            compound = json.loads(line)
            table[compound['name']] = compound['id']
            if 'synonyms' in compound:
                for synonym in compound['synonyms']:
                    table[synonym] = compound['id']
    return table

if __name__ == "__main__":
    table = parse_table(sys.argv[1])
    print(table)
