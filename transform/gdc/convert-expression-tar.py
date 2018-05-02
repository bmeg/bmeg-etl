#!/usr/bin/env python

import os
import re
import sys
import csv
import gzip
import json
import string
import argparse
import tempfile
import subprocess
from sets import Set
from pprint import pprint
from google.protobuf import json_format
import rna_pb2
from google.protobuf import json_format

import tarfile

def process_expression(data, source, project_id, sample_id, emit):
    out = rna_pb2.GeneExpression()
    for line in data:
        row = line.rstrip().split("\t")
        k = row[0]
        v = float(row[1])
        if v != 0.0:
            out.expressions[k] = v
    out.source = source
    out.id = sample_id
    out.biosample_id = sample_id
    emit(out)

def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    return json.dumps(msg)

def convert(args):
    sample_map = {}
    project_map = {}
    if args.filemap is not None:
        with open(args.filemap) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample_map[row['file_id']] = row["sample_id"]
                project_map[row["file_id"]] = row["project_id"]
    
    handle = open(args.out, "w")
    def emit(record):
        handle.write(message_to_json(record) + "\n")

    tmpdir = tempfile.mkdtemp(dir=".", prefix="gdc_extract_")
    subprocess.check_call("tar xvzf %s" % (os.path.abspath(args.tar)), shell=True, cwd=tmpdir)

    with open(os.path.join(tmpdir, "MANIFEST.txt")) as man:
        file_map = {}
        for line in man:
            row = line.rstrip().split("\t")
            file_map[row[0]] = row[1]
    
    for k, v in file_map.items():
        p = os.path.join(tmpdir, v)
        if os.path.exists(p):
            with gzip.GzipFile(p) as f:
                project = ""
                sample = k
                if k in sample_map:
                    sample = sample_map[k]
                if k in project_map:
                    project = project_map[k]
                process_expression(f, args.source, project, sample, emit)

    handle.close()

def parse_args(args):
    args = args[1:]
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', type=str, default="test.out", help='path to output file')
    parser.add_argument("--filemap", default=None)
    parser.add_argument("--source", default='gdc')
    parser.add_argument("tar")
    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert(options)
