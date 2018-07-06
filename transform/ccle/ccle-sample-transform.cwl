
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest

baseCommand:
  - /opt/convert-ccle-sample.py

arguments:
  - valueFrom: "Sample.json"
    position: 2

inputs:
  SAMPLE_TSV:
    type: File
    inputBinding:
      position: 1

outputs:
  OUTPUT:
    type: File
    outputBinding:
      glob: "Sample.json"
