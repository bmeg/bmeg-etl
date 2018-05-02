#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/gdc-transform:latest


baseCommand: 
  - python
  - /opt/convert-expression-tar.py

arguments:
  - "--out"
  - $(inputs.outfile)

inputs:
  outfile:
    type: string
    default: rna_seq.json

  FILE_MAP:
    type: File
    inputBinding:
      prefix: "--filemap"
  TARBALL:
    type: File
    inputBinding:
      position: 2

outputs:
  DATA:
    type: File
    outputBinding:
      glob: $(inputs.outfile)
