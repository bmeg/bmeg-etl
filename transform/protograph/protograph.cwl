
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: spanglry/protograph
baseCommand:
  - "java"
  - "-jar"
  - "/protograph.jar"
  - "--output"
  - "out"

inputs:
  PROTOGRAPH:
    type: File
    inputBinding:
      prefix: "--protograph"
  LABEL:
    type: string
    inputBinding:
      prefix: "--label"
  INPUT:
    type: File
    inputBinding:
      prefix: "--input"

outputs:
  VERTEXES:
    type: File
    outputBinding:
      glob: out.Vertex.json
  EDGES:
    type: File
    outputBinding:
      glob: out.Edge.json
