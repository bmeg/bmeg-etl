#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/ctdd-transform:latest


baseCommand: 
  - python
  - /opt/convert-ctdd.py
  - "--format"
  - json
  - "--multi"
  - ctdd
  - "--pubchem"
  - "/opt/ctdd_pubchem.table"

inputs:
  CTRP_ZIP:
    type: File
    inputBinding: 
      prefix: "--zip"

outputs:
  RESPONSE:
    type: File
    outputBinding:
      glob: ctdd.bmeg.ResponseCurve.json
