#!/usr/bin/env python

import json
import sys
import re

re_section = re.compile(r'^\[(.*)\]')
re_field = re.compile(r'^(\w+): (.*)$')

def obo_parse(handle):
    rec = None
    for line in handle:
        res = re_section.search(line)
        if res:
            if rec is not None:
                yield rec
            rec = None
            if res.group(1) == "Term":
                rec = {"type": res.group(1)}
        else:
            if rec is not None:
                res = re_field.search(line)
                if res:
                    key = res.group(1)
                    val = res.group(2)
                    val = val.split(" ! ")[0]
                    if key in rec:
                        rec[key].append(val)
                    else:
                        rec[key] = [val]
    if rec is not None:
        yield rec

def transform(obo_file, output_base):
    with open(output_base, "wt") as ohandle:
        with open(obo_file) as handle:
            for rec in obo_parse(handle):
                ohandle.write(json.dumps(rec))
                ohandle.write("\n")

if __name__ == "__main__":
    transform(sys.argv[1], sys.argv[2])
