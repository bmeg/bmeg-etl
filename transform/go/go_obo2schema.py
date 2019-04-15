#!/usr/bin/env python

import re
import sys

from bmeg import GeneOntologyTerm
from bmeg.emitter import JSONEmitter

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


def unquote(s):
    res = re.search(r'"(.*)"', s)
    if res:
        return res.group(1)
    return s


if __name__ == "__main__":

    emitter = JSONEmitter(sys.argv[2])

    with open(sys.argv[1]) as handle:
        for rec in obo_parse(handle):
            go_id = rec['id'][0]
            go_name = rec['name'][0]
            go_namespace = rec['namespace'][0]
            go_definition = unquote(rec['def'][0])
            synonym = []
            is_a = []
            xref = []
            if 'synonym' in rec:
                for i in rec['synonym']:
                    synonym.append(unquote(i))
            if 'is_a' in rec:
                for i in rec['is_a']:
                    is_a.append(i)
            if 'xref' in rec:
                for i in rec['xref']:
                    xref.append(i.split(" ")[0])
            emitter.emit_vertex(GeneOntologyTerm(
                go_id=go_id, name=go_name,
                definition=go_definition, namespace=go_namespace,
                synonym=synonym, xref=xref, is_a=is_a
            ))

    emitter.close()
