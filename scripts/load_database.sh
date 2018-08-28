#!/bin/bash

graph="bmeg-test"
arachne drop $graph
arachne create $graph

gofast="--numInsertionWorkers 8 --writeConcern 0 --bypassDocumentValidation"

for f in $(ls -1 outputs/*/*.Vertex.json); do
  mongoimport -d arachnedb -c ${graph}_vertices --type json --file $f $gofast
done

for f in $(ls -1 outputs/*/*.Edge.json); do
  mongoimport -d arachnedb -c ${graph}_edges --type json --file $f $gofast
done
