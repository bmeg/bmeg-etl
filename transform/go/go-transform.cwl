

cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest


baseCommand:
  - python
  - /opt/go_obo2schema.py

inputs:
  OBO:
    type: File
    inputBinding:
      position: 1

outputs:
  GO_JSON:
    type: stdout

stdout: go.json
