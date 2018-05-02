
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/pubchem-transform:latest
baseCommand:
  - /opt/pubchem-extract.py
  - compound-extract

inputs:
  OUTNAME:
    type: string
    default: compound.json
    inputBinding:
      position: 1
      prefix: "--out"

  CIDS:
    type: string[]
    inputBinding:
      position: 2

outputs:
  CASE_LIST:
    type: File
    outputBinding:
      glob: $(inputs.OUTNAME)
