
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/ccle-transform:latest

baseCommand:
  - /command/convert-ccle-sample.py
  
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
