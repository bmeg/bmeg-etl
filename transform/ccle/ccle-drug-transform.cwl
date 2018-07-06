
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest

baseCommand:
  - /opt/convert-ccle.py
  - "--format"
  - json
  - "--pubchem"
  - /in/ccle_pubchem.txt
  - "--multi"
  - "ccle"

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
    type: File
    outputBinding:
      glob: ccle.bmeg.ResponseCurve.json
