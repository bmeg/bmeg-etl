


cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: appropriate/curl

baseCommand:
  - curl

inputs:
  URL:
    type: string
    inputBinding:
      position: 2
  NAME:
    type: string
    inputBinding:
      position: 1
      prefix: "-o"

outputs:
  OUTPUT:
    type: File
    outputBinding:
      glob: $(inputs.NAME)
      