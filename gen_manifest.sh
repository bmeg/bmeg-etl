#!/bin/bash

find ./outputs -name "*.json.gz" | grep -v Compound | grep -v Allele > bmeg_file_manifest.txt 
find ./outputs/compound -name "normalized.*.json.gz" >> bmeg_file_manifest.txt
find ./outputs/allele -name "normalized.*.json.gz" >> bmeg_file_manifest.txt
