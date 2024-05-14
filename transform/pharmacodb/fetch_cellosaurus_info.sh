#!/bin/bash


fetch_cell_info() {
	# name is sample_id 
	# accession_id is cellosaurus id
	# diseases is cellosaurus disease
	local output_file="cell_lines_cellosaurus_mapping.json"
	    curl -s -X POST https://pharmacodb.ca/graphql \
        -H "Content-Type: application/json" \
        -d '{"query":"{ cell_line(: '"$cellId"') { id name accession_id diseases tissue { id name } synonyms { name  dataset { id name } } } "}' \
        -o "$output_file"
}
