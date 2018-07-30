#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys
import tempfile
import yaml


def run_run(args):
    args.exec_dir = os.path.abspath(args.exec_dir)

    with open(args.cwl) as handle:
        workflow = yaml.load(handle.read())

    if 'bmeg:split' in workflow:
        # if they have defined an input file as a split, create an array of
        # existing variable inputs
        inputArray = []
        splitFile = os.path.abspath(os.path.join(
            args.build_dir, workflow["bmeg:split"]))
        with open(splitFile) as handle:
            for line in handle:
                inputArray.append(json.loads(line))
    else:
        # otherwise start a single run, with no input varibles
        inputArray = [{}]

    for buildInputs in inputArray:
        inputs = buildInputs
        if len(workflow['inputs']):
            for k, v in workflow['inputs'].items():
                if 'bmeg:key' in v:
                    inputPath = os.path.abspath(os.path.join(
                        args.build_dir, v['bmeg:key'].format(**buildInputs)))
                    inputs[k] = {"class": "File", "path": inputPath}

        outputNotFound = False
        outMapping = {}
        for k, v in workflow['outputs'].items():
            if 'bmeg:key' in v:
                path = os.path.abspath(os.path.join(
                    args.build_dir, v['bmeg:key'].format(**buildInputs)))
                outMapping[k] = path
                if not os.path.exists(path):
                    outputNotFound = True

        if outputNotFound:

            if not os.path.exists(os.path.join(args.exec_dir, "inputs")):
                os.makedirs(os.path.join(args.exec_dir, "inputs"))

            cwl = os.path.splitext(os.path.basename(args.cwl))[0]
            inputFile = tempfile.NamedTemporaryFile(
                dir=os.path.join(args.exec_dir, "inputs"),
                prefix="bmegbuild-" + cwl + "-",
                suffix=".json",
                delete=False
            )
            inputFile.write(json.dumps(inputs))
            inputFile.close()

            tmpdir = os.path.join(args.exec_dir, "tmp/tmp")
            outdir = os.path.join(args.exec_dir, "out")
            cachedir = os.path.join(args.exec_dir, "cache")

            try:
                cmd = [
                    "cwltool",
                    "--outdir", outdir,
                    "--tmpdir-prefix", tmpdir,
                    "--cachedir", cachedir,
                    args.cwl, inputFile.name
                ]
                print("Running:", " ".join(cmd))
                print("Inputs:", inputs)
                if args.dry_run:
                    continue
                outJSON = subprocess.check_output(cmd)
                outData = json.loads(outJSON)

                for k, v in outData.items():
                    if k in outMapping:
                        outPath = outMapping[k]
                        if not os.path.exists(os.path.dirname(outPath)):
                            os.makedirs(os.path.dirname(outPath))
                        print("Moving", v['path'], outPath)
                        os.rename(v['path'], outPath)
            except Exception as e:
                raise e
            finally:
                os.remove(inputFile.name)

        else:
            print("All outputs found")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd", help='sub-command help')

    parser_run = subparsers.add_parser('run', help='run a BMEG CWL ETL tool')
    parser_run.set_defaults(func=run_run)
    parser_run.add_argument("cwl", help='CWL file')
    parser_run.add_argument("build_dir", help='BMEG build directory')
    parser_run.add_argument("--dry-run", action="store_true",
                            help="print generated CWL inputs and exit")
    parser_run.add_argument("--exec-dir", default=".",
                            help="base directory to use for cwltool '--outdir', \
                            '--tmpdir-prefix', '--tmp-outdir-prefix', \
                            '--cachedir'")
    args = parser.parse_args()
    if not args.cmd:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args.func(args)
