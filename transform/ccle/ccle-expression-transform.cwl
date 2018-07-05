
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest

baseCommand:
  - /opt/convert-ccle-expression.py

inputs:
  GCT:
    type: File
    inputBinding:
      position: 1

outputs:
  OUTPUT:
    type: File
    outputBinding:
      glob: expression_ccle.json
