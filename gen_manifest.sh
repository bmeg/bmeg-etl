#!/bin/bash

find ./outputs -name "*.json.gz" | grep -v Compound > bmeg_file_manifest.txt 
find ./outputs/compound -name "normalized.*.json.gz" >> bmeg_file_manifest.txt
