#!/usr/bin/env python

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

    subprocess.check_call(["cwltool", args.cwl, inputFile.name])



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
