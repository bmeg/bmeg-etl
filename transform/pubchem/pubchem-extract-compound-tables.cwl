
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/pubchem-transform:latest
baseCommand:
  - /opt/pubchem-extract.py
  - table-xform

inputs:
  OUTNAME:
    type: string
    default: compound-tables.json
    inputBinding:
      position: 1
      prefix: "--table_xform_out"


outputs:
  CASE_LIST:
    type: File
    outputBinding:
      glob: $(inputs.OUTNAME)
