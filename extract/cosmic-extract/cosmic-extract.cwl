#!/usr/bin/env cwl-runner

# this should produce a docker run of the form:
# docker run --rm -it -v /tmp/cosmic:/output -v ~/g2p-aggregator/harvester/CosmicMutantExport.tsv.gz:/g2p-aggregator/harvester/CosmicMutantExport.tsv.gz  cosmic-extract
# note: the swiftPath attribute is somewhat ad-hoc, there is no official support for openstack swift

cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: cosmic-extract:latest
baseCommand: /g2p-aggregator/harvester/make-all.sh
inputs:
  COSMIC:
    type: File
    inputBinding:
      position: 1
outputs:
  output:
    type: Directory
    outputBinding:
      glob: output