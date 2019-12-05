#!/bin/bash

set -e

if [ "$#" -ne 2 ]; then
		printf "Illegal number of parameters.\n\n"
		printf "Usage: rsync_manifest_files.sh <file_manifest> [USER@]HOST:DEST\n"
		printf "       rsync_manifest_files.sh <file_manifest> DEST\n"
		exit 1
fi

file_manifest=$(realpath $1)
dst=$2

echo $PWD
for f in $(cat $file_manifest); do
		if [ ! -f $f ]; then
				echo "file $f does not exist!"
		else
				rsync $f ${dst}/$(dirname $f)/
		fi
done
