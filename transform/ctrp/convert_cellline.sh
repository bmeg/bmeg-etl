#!/bin/bash

curl -o DepMap-2018q4-celllines.csv -L https://depmap.org/portal/download/api/download/external?file_name=processed_portal_downloads%2Fdepmap-public-cell-line-metadata-183e.1%2FDepMap-2018q4-celllines.csv

cat DepMap-2018q4-celllines.csv | grep -v "NORMAL" | grep -v MERGED_TO_ACH | grep -v "TT_" | \
    grep -v MS1_ | grep -v KMH2_THYROID | tr "," "_" | awk -F "_" '{print $2 "\t" $1}' > ctrp-cellline.table
