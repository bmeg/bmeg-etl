#!/usr/bin/env python

import os
import sys
import csv
import json
import argparse
import tempfile
import subprocess
from google.protobuf import json_format
from bmeg.methylation_pb2 import MethylationPanel

def meth_parse_new(path):
    with open(path) as handle:
        out = {}
        reader = csv.DictReader(handle, delimiter="\t")
        out = MethylationPanel()
        for row in reader:
            #print row
            try:
                out.beta_value[row['Composite Element REF']] = float(row['Beta_value'])
            except ValueError:
                pass
        return out


def message_to_json(message):
    msg = json_format.MessageToDict(message)
    return json.dumps(msg)

def parse_tar(args):

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
    subprocess.check_call("tar xzf %s" % (os.path.abspath(args.tar)), shell=True, cwd=tmpdir)

    with open(os.path.join(tmpdir, "MANIFEST.txt")) as man:
        file_map = {}
        for line in man:
            row = line.rstrip().split("\t")
            file_map[row[0]] = os.path.join(tmpdir, row[1])

    for k, v in file_map.items():
        if os.path.exists(v):
            o = meth_parse_new(v)
            o.biosample_id = sample_map[k]
            emit(o)

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--file", default=None)
  parser.add_argument("--tar", default=None)
  parser.add_argument("--filemap", default=None)
  parser.add_argument("--source", default='gdc')
  parser.add_argument("--out", default="methylation.data.json")
  args = parser.parse_args()

  if args.tar is not None:
      parse_tar(args)
