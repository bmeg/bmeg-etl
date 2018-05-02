
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/pubmed-transform:latest

baseCommand:
  - python
  - /opt/pubmed.py

inputs:
  file:
    type: File
    inputBinding:
      position: 1

stdout: pubmed.json

outputs:
  DATA:
    type: stdout
