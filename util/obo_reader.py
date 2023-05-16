#!/usr/bin/env python

import re
import sys
import json

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
                    val = re.split(" ! | \(|\)", val)
                    val = ":".join(val[0:3])
                    if key in rec:
                        rec[key].append(val)
                    else:
                        rec[key] = [val]

    if rec is not None:
        yield rec


def unquote(s):
    res = re.search(r'"(.*)"', s)
    if res:
        return res.group(1)
    return s


with open(sys.argv[1]) as handle:
    for rec in obo_parse(handle):
        print(json.dumps(rec))
