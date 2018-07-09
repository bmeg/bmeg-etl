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

    if 'bmeg:split' in workflow:
        #if they have defined an input file as a split, create an array of existing
        #variable inputs
        inputArray = []
        with open(workflow["bmeg:split"]) as handle:
            for line in handle:
                inputArray.append(json.loads(line))
    else:
        #otherwise start a single run, with no input varibles
        inputArray = [{}]

    for buildInputs in inputArray:
        inputs = buildInputs
        if len(workflow['inputs']):
            for k,v in workflow['inputs'].items():
                if 'bmeg:key' in v:
                    inputPath = v['bmeg:key'].format(**buildInputs)
                    inputs[k] = {"class":"File", "path" : inputPath}

        outputNotFound = False
        outMapping = {}
        for k,v in workflow['outputs'].items():
            if 'bmeg:key' in v:
                path = v['bmeg:key'].format(**buildInputs)
                outMapping[k] = path
                if not os.path.exists(path):
                    outputNotFound = True

        if outputNotFound:
            inputFile = tempfile.NamedTemporaryFile(dir="./", prefix="bmegbuild-", suffix=".json", delete=False)
            inputFile.write(json.dumps(inputs))
            inputFile.close()
            print inputFile.name

            outJSON = subprocess.check_output(["cwltool", args.cwl, inputFile.name])
            outData = json.loads(outJSON)


            for k, v in outData.items():
                if k in outMapping:
                    outPath = outMapping[k]
                    if not os.path.exists(os.path.dirname(outPath)):
                        os.makedirs(os.path.dirname(outPath))
                    print("Moving", v['path'], outPath)
                    os.rename(v['path'], outPath)
            os.unlink(inputFile.name)
        else:
            print("All outputs found")



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
