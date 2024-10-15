#!/bin/bash

fetch_experiments() {
    local page=$1
    local per_page=$2
    local output_file=$3
    curl -s -X POST https://pharmacodb.ca/graphql \
        -H "Content-Type: application/json" \
        -d '{"query":"{ experiments(page: '"$page"', per_page: '"$per_page"') { id cell_line { id name tissue { id name } } compound { id name annotation { smiles inchikey pubchem chembl fda_status } } tissue { id name } dataset { id name } profile { HS Einf EC50 AAC IC50 DSS1 DSS2 DSS3 } dose_response { dose response } } }"}' | \
        jq -c '.data.experiments'  >> "$output_file"
    echo "fetched page:  $page"
}

fetch_all_experiments() {
    local per_page=10000
    local output_file="../../source/pharmacodb/graphql_dump/experiments.ndjson"
    # range via: https://pharmacodb.ca/graphql 
    for ((page=1; page<=632; page++))
    do
        fetch_experiments $page $per_page $output_file
    done
}

fetch_all_experimen