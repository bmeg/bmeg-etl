
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: pfam-transform:latest

baseCommand:
  - python
  - /opt/pfam_transform.py

arguments:
  - "file"

inputs:
  files:
    type: File[]
    inputBinding:
      position: 2

stdout: pfam.json

outputs:
  DATA:
    type: stdout