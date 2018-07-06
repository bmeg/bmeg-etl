
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/gdc-extract:latest
baseCommand:
  - /opt/gdc-scan.py
  - "--out"
  - gdc_case_scan.json
  - cases
  - list

inputs:
  blank:
    type: [boolean, "null"]
outputs:
  CASE_LIST:
    type: File
    outputBinding:
      glob: gdc_case_scan.json
