#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/gdc-transform:latest

baseCommand:
  - python
  - /opt/convert-clinical.py

arguments:
  - "--out"
  - ./

inputs:
  TARBALL:
    type: File
    inputBinding:
      prefix: "--tar"

outputs:
  BIOSAMPLE:
    type: File
    outputBinding:
      glob: tcga.Biosample.json
