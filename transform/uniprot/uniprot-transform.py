#!/usr/bin/env python3

import argparse
import gzip
import json
import logging
import re
import sys
import xml.sax
import os
from ftplib import FTP

from bmeg.vertex import Protein
from bmeg.emitter import JSONEmitter



class StackLevel:
    def __init__(self, name, attrs):
        self.data = {}
        self.name = name
        self.has_children = False
        self.attrs = attrs

def string_pass(e, v, attrs):
    return v

def create_list(e, v, attrs):
    return [v]

def emit_entry(e, v, attrs, **kwds):
    #print(kwds)
    pass

def pass_attrs(e, v, attrs, **kwds):
    attrs.update(kwds)
    return attrs

def create_dict_list(e, v, attrs, **kwds):
    return [kwds]

def dbRefExtract(e, v, attrs, **kwds):
    print(attrs, kwds)

f_map = [
    (["uniprot","entry","accession"], None, create_list),
    (["uniprot","entry"], None, emit_entry),
    (["uniprot","entry","dbReference"], None, dbRefExtract),
    (["uniprot","entry","dbReference", "property"], None, create_dict_list),
]

NOT_FOUND = set()

class UniprotHandler(xml.sax.ContentHandler):
    def __init__(self, record_write):
        xml.sax.ContentHandler.__init__(self)
        self.record_write = record_write
        self.stack = []

    def startElement(self, name, attrs):
        if len(self.stack):
            self.stack[-1].has_children = True
        self.stack.append(StackLevel(name, dict(attrs.items())))
        self.buffer = ""

    def characters(self, text):
        self.buffer += text

    def endElement(self, name):
        stack_id = list(i.name for i in self.stack)
        level = self.stack.pop()
        found = False
        for s, out_name, f in f_map:
            if stack_match(s, stack_id):
                if out_name is None:
                    out_name = stack_id[-1]
                v = f(self.record_write, self.buffer, level.attrs, **level.data)
                if v is not None:
                    if isinstance(v, list):
                        if out_name in self.stack[-1].data:
                            self.stack[-1].data[out_name].extend(v)
                        else:
                            self.stack[-1].data[out_name] = v
                    elif isinstance(v, dict):
                        if out_name in self.stack[-1].data:
                            self.stack[-1].data[out_name] = dict(self.stack[-1].data, **v)
                        else:
                            self.stack[-1].data[out_name] = v
                    else:
                        self.stack[-1].data[out_name] = v
        if not found:
            n = ",".join(stack_id)
            if n not in NOT_FOUND:
                NOT_FOUND.add(n)
                logging.warning("combiner for %s not found" % (",".join(stack_id)))
        self.buffer = ""


def stack_match(query, elem):
    if len(query) != len(elem):
        return False
    for q, e in zip(query, elem):
        if isinstance(q, list):
            if e not in q:
                return False
        else:
            if q != "*" and q != e:
                return False
    return True


def parse_uniprot(handle, emitter):
    handler = UniprotHandler(emitter)
    parser = xml.sax.make_parser()
    parser.setContentHandler(handler)
    parser.parse(handle)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    parser.add_argument("--output", default="uniprot")
    args = parser.parse_args()

    emitter = JSONEmitter(args.output)
    with open(args.file) as handle:
        parse_uniprot(handle, emitter)
