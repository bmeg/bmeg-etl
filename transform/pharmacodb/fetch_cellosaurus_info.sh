#!/bin/bash

fetch_cell_info() {
    local celline_id=$1
    echo "in fetch_cell_info - fetching info for cell-line ID: $celline_id"
    curl -s -X POST https://pharmacodb.ca/graphql \
        -H "Content-Type: application/json" \
        -d '{"query":"{ cell_line(cellId: '"$celline_id"') { id name accession_id diseases tissue { id name } synonyms { name dataset { id name } } } }"}' | \
        jq -c '.data.cell_line' >> "../../source/pharmacodb/graphql_dump/cells_cellosaurus_info.ndjson"
}

fetch_all_cell_line_info() {
    > "../../source/pharmacodb/graphql_dump/cells_cellosaurus_info.ndjson"

    jq -r '.data.cell_lines[].id' ../../source/pharmacodb/graphql_dump/cell_lines.json | while read celline_id; do
        echo "iterating over cell-line ID: $celline_id"
        fetch_cell_info $celline_id
    done
}

fetch_all_cell_line_info
