

cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest

baseCommand:
  - "python3.7"
  - "/opt/transform-rna-matrix"
  - "--output-prefix"
  - "rna"

inputs:
  matrix:
    type: File
    inputBinding:
      prefix: "--matrix"

  gzip:
    type: boolean
    inputBinding:
      prefix: "--gz"

outputs:
  GeneExpression:
    type: File
    outputBinding:
      glob: rna.GeneExpression.Vertex.json
  ExpressionOf:
    type: File
    outputBinding:
      glob: rna.ExpressionOf.Edge.json
