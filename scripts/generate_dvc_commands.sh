#!/bin/bash

set -e 

dvc pipeline show outputs.bmeg_manifest.dvc | sort -r | python ./transform/dvc/dvc2cmd.py - > dvc_commands.txt
