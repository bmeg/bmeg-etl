
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: gdc-extract:latest
baseCommand:
  - /opt/gdc-scan.py
  - files
  - download

inputs:
  ID:
    type: string
    inputBinding:
      prefix: "--id"
outputs:
  FILE:
    type: File
    outputBinding:
      glob: "*.out"

