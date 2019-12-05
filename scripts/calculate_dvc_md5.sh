#!/bin/bash

set -e 

yaml2json() {
		python -c 'import sys, yaml, json, dvc.utils; d = yaml.load(sys.stdin, Loader=yaml.FullLoader); print(dvc.utils.dict_md5(d))'
}
cat $1 | grep -v "^md5" | grep -v "persist" | grep -v "metric" | grep -v "tags" | yaml2json
