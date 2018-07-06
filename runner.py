#!/usr/bin/env python

import os
import yaml
import json
import argparse
import tempfile
import subprocess

def run_list(args):
    print "Hello"

def run_run(args):
    with open(args.cwl) as handle:
        workflow = yaml.load(handle.read())

    inputs = {}
    if len(workflow['inputs']):
        for k,v in workflow['inputs'].items():
            if 'bmeg:key' in v:
                inputs[k] = {"class":"File", "path" : v['bmeg:key']}

    inputFile = tempfile.NamedTemporaryFile(dir="./", prefix="bmegbuild-", suffix=".json", delete=False)
    inputFile.write(json.dumps(inputs))
    inputFile.close()
    print inputFile.name

    outJSON = subprocess.check_output(["cwltool", args.cwl, inputFile.name])
    outData = json.loads(outJSON)

    outMapping = {}
    for k,v in workflow['outputs'].items():
        if 'bmeg:key' in v:
            outMapping[k] = v['bmeg:key']

    for k, v in outData.items():
        if k in outMapping:
            outPath = outMapping[k]
            if not os.path.exists(os.path.dirname(outPath)):
                os.makedirs(os.path.dirname(outPath))
            os.rename(v['path'], outPath)
    os.unlink(inputFile)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')
    parser_list = subparsers.add_parser('list', help='a help')
    parser_list.set_defaults(func=run_list)

    parser_run = subparsers.add_parser('run', help='a help')
    parser_run.set_defaults(func=run_run)
    parser_run.add_argument("cwl")

    args = parser.parse_args()

    args.func(args)
