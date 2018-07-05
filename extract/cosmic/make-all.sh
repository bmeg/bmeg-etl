#!/bin/bash

ln -s $1 ./CosmicMutantExport.tsv.gz
cp -r /g2p-aggregator/harvester/* ./
make SHELL=/bin/bash
mkdir output
mv *.tsv output/