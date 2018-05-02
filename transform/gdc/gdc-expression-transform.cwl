#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: gdc-transform:latest


baseCommand: 
  - python
  - /opt/convert-expression.py
  - "--gz"

arguments:
  - "--out"
  - $(inputs.outfile)

inputs:
  outfile:
    type: string
    default: rna_seq.json

  PROJECT:
    type: string
    inputBinding:
      prefix: "--project"
  SAMPLE:
    type: string
    inputBinding:
      prefix: "--sample"
      position: 2
  FILE:
    type: File
    inputBinding:
      position: 3

outputs:
  DATA:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.outfile)
