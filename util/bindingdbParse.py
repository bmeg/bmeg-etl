#!/usr/bin/env python

import re
import sys
import json

#This parser is designed specifically for parsing bindingDB data. It is built around the fact that bindingDB repeats the last 12 columns a certain number of times based on how many peptides are in the complex described. If bindingDB changes the number of columnns in their table, a change to the indexes used to loop woill likely be needed

def tsv_parse(handle):
    #pruning off repetitive columns
    columns = handle.readline().strip('\n').split("\t")[:49]

    for line in handle:
        values = line.strip('\n').split("\t")
        #When there are multiple protiens, there are multiple columns. Combining identical columns into individual lists for each field in each row
        for i in range(0,12):
            values[37+i] = values[37+i::12]
        values = values[:49]
        #making sure it's a normal object before passing it on
        if len(values) == 49:
            row = dict(zip(columns,values))
            yield row

with open(sys.argv[1]) as handle:
    for rec in tsv_parse(handle):
        print(json.dumps(rec))
