#!/bin/bash

wget -O source/pathway_commons/pathways.txt.gz https://www.pathwaycommons.org/archives/PC2/v11/pathways.txt.gz
wget -O source/pathway_commons/paxtools.jar https://www.pathwaycommons.org/archives/PC2/v11/paxtools.jar
wget -O source/pathway_commons/PathwayCommons11.Detailed.BIOPAX.owl.gz https://www.pathwaycommons.org/archives/PC2/v11/PathwayCommons11.All.BIOPAX.owl.gz
java -Xmx50g -jar source/pathway_commons/paxtools.jar toSIF source/pathway_commons/PathwayCommons11.All.BIOPAX.owl.gz source/pathway_commons/pc11.all.sif seqDb=hgnc exclude=neighbor_of MEDIATOR PUBMED PMC COMMENTS PATHWAY PATHWAY_URI RESOURCE SOURCE_LOC TARGET_LOC
