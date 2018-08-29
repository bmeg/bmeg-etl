
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest

baseCommand:
  - python3.7
  - /opt/transform/pubmed/pubmed.py

inputs:
  file:
    type: File
    inputBinding:
      position: 1

outputs:
  DATA:
    type: File
    outputBinding:
      glob: outputs/pubmed/pubmed.Pubmed.Vertex.json
