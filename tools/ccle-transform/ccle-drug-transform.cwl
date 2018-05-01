
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/ccle-transform:latest

baseCommand:
  - /command/convert-ccle.py
  - "--format"
  - json
  - "--pubchem"
  - /in/ccle_pubchem.txt


arguments:
  - "--multi"
  - "$(inputs.OUTPATH)"

inputs:

  DRUG:
    type: File
    inputBinding:
      prefix: "--drug"
  OUTPATH:
    type: [string, "null"]
    default: "ccle-data"

outputs:
  OUTPUT:
    type: File[]
    outputBinding:
      glob: "*.json"
