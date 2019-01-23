#!/bin/bash

set -e

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ "$#" -ne 1 ]; then
		printf "Usage:\n	generate_bmeg_file_manifest.sh <directory to search>\n\n"
		exit 1
fi

DIR=$1

VERTEX_FILES=($DIR/outputs/**/*.Vertex.json.gz)
EDGE_FILES=($DIR/outputs/**/*.Edge.json.gz)
FILES=(${VERTEX_FILES[@]} ${EDGE_FILES[@]})

EXCEPTIONS=()

echo "generating DVC command..."

DVC_CMD="dvc run --file outputs.bmeg_manifest.dvc --yes "

for f in ${FILES[@]}; do
		if [[ ! " ${EXCEPTIONS[@]} " =~ " ${f} " ]]; then
				DVC_CMD+=" -d $f "
		else
				echo "excluding $f..."
		fi
done

DVC_CMD+=" \"python3 transform/dvc/bmeg_file_manifest.py\""

echo "running DVC..."

# echo $DVC_CMD | bash
