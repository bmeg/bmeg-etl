#!/bin/bash

set -e

VERTEX_FILES=(outputs/**/*.Vertex.json.gz)
EDGE_FILES=(outputs/**/*.Edge.json.gz)
FILES=(${VERTEX_FILES[@]} ${EDGE_FILES[@]})

EXCEPTIONS=(
		"outputs/ccle/maf.Allele.Vertex.json.gz"
		"outputs/ccle/drug_response.Compound.Vertex.json.gz"
		"outputs/ccle/drug_response.ResponseTo.Edge.json.gz"
		"outputs/ctrp/ResponseTo.Edge.json.gz"
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
)

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
echo $DVC_CMD
echo $DVC_CMD | bash
