
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: pfam-transform:latest

baseCommand:
  - python
  - /opt/pfam_transform.py


arguments:
  - "download"

inputs:
  all:
    type: boolean?
    default: False
    inputBinding:
      prefix: "--all"

  ids:
    type: ["null", "string[]"]
    inputBinding:
      position: 2

outputs:
  DATA:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.xml"
