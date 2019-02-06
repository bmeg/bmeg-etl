#!/bin/bash

set -e

DVC_FILES=( $(find . -maxdepth 1 -type f -name 'outputs.*.dvc' | sort) )
FILES=()
for FILE in ${DVC_FILES[@]}; do
		if [[ ! "./outputs.bmeg_manifest.dvc" == "${FILE}" ]]; then
				OUT=($(cat $FILE | grep "path:" | grep -v "source/" | sed 's~  path: ~~g' | sort))
				for F in ${OUT[@]}; do
						if [ -d "${F}" ]; then
								VERTEX_FILES=( $(find ${F} -type f -name "*.Vertex.json.gz" | sort) )
								EDGE_FILES=( $(find ${F} -type f -name "*.Edge.json.gz" | sort) )
								FILES=(${FILES[@]} ${VERTEX_FILES[@]} ${EDGE_FILES[@]})
						else
								FILES=(${FILES[@]} $F)
						fi
				done
		fi
done

EXCEPTIONS=(
		"outputs/ccle/maf.Allele.Vertex.json.gz"
		"outputs/ccle/drug_response.Compound.Vertex.json.gz"
		"outputs/ccle/drug_response.ResponseTo.Edge.json.gz"
		"outputs/ctrp/ctrp.Compound.Vertex.json.gz"
		"outputs/ctrp/ctrp.ResponseTo.Edge.json.gz"
		"outputs/g2p/Allele.Vertex.json.gz"
		"outputs/g2p/Compound.Vertex.json.gz"
		"outputs/g2p/HasEnvironment.Edge.json.gz"
		"outputs/gdc/Compound.Vertex.json.gz"
		"outputs/gdc/TreatedWith.Edge.json.gz"
		"outputs/gdsc/gdsc.Compound.Vertex.json.gz"
		"outputs/gdsc/gdsc.ResponseTo.Edge.json.gz"
		"outputs/mc3/Allele.Vertex.json.gz"
		"outputs/meta/Command.Vertex.json.gz"
		"outputs/meta/File.Vertex.json.gz"
		"outputs/meta/Reads.Edge.json.gz"
		"outputs/meta/Writes.Edge.json.gz"
		"outputs/meta/bmeg_file_manifest.txt"
)

echo "generating DVC command..."

DVC_CMD="dvc run --file outputs.bmeg_manifest.dvc --yes --ignore-build-cache"

for FILE in ${FILES[@]}; do
		if [[ ! " ${EXCEPTIONS[@]} " =~ " ${FILE} " ]]; then
				DVC_CMD+=" -d $FILE "
		else
				echo "excluding $FILE..."
		fi
done

DVC_CMD+=" \"python3 transform/dvc/bmeg_file_manifest.py\""

echo "running DVC..."
echo $DVC_CMD
echo $DVC_CMD | bash
