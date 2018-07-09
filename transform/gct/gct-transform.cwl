
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest

baseCommand:
  - "/opt/gct"

inputs:
  GCT:
    type: File
    inputBinding:
      prefix: "-filepath"

  SOURCE:
    type: string
    inputBinding:
      prefix: "-source"

  SCALE:
    type: string
    inputBinding:
      prefix: "-scale"

  GZIPPED:
    type: boolean?
    inputBinding:
      prefix: "-gzipped"

outputs:
  OUTPUT:
    type: File
    outputBinding:
      glob: out.json
