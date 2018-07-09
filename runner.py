#!/usr/bin/env python

from __future__ import print_function

import os
import yaml
import json
import argparse
import tempfile
import subprocess


def run_run(args):
    with open(args.cwl) as handle:
        workflow = yaml.load(handle.read())

    inputs = {}
    if len(workflow['inputs']):
        for k,v in workflow['inputs'].items():
            if 'bmeg:key' in v:
                inputs[k] = {"class":"File", "path" : os.path.join(args.build_dir, v['bmeg:key'])}

    if args.dry_run:
        print(inputs)
        return

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

    tmpdir = os.path.join(args.exec_dir, "tmp")
    outdir = os.path.join(args.exec_dir, "out")
    cachedir = os.path.join(args.exec_dir, "cache")

    try:
        outJSON = subprocess.check_output([
            "cwltool",
            "--outdir", outdir,
            "--tmpdir-prefix", tmpdir,
            "--cachedir", cachedir,
            args.cwl, inputFile.name
        ])
        outData = json.loads(outJSON)

        outMapping = {}
        for k, v in workflow['outputs'].items():
            if 'bmeg:key' in v:
                outMapping[k] = os.path.join(args.build_dir, v['bmeg:key'])

        for k, v in outData.items():
            if k in outMapping:
                outPath = outMapping[k]
                if not os.path.exists(os.path.dirname(outPath)):
                    os.makedirs(os.path.dirname(outPath))
                os.rename(v['path'], outPath)
        print(outData)

    except Exception as e:
        raise e

    finally:
        os.remove(inputFile.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')
    
    parser_run = subparsers.add_parser('run', help='run a BMEG CWL ETL tool')
    parser_run.set_defaults(func=run_run)
    parser_run.add_argument("cwl", help='CWL file')
    parser_run.add_argument("build_dir", help='BMEG build directory')
    parser_run.add_argument("--dry-run", action="store_true", help="print generated CWL inputs and exit")
    parser_run.add_argument("--exec-dir", default=".", help="base directory to use for cwltool '--outdir', '--tmpdir-prefix', '--tmp-outdir-prefix', '--cachedir'")

    args = parser.parse_args()
    args.exec_dir = os.path.abspath(args.exec_dir)
    args.func(args)
