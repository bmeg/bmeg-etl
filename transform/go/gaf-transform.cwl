

cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest


baseCommand:
  - python
  - /opt/go_gaf2schema.py

arguments:
  - valueFrom: gaf.json
    position: 3

inputs:
  GAF:
    type: File
    inputBinding:
      position: 1
  UNIPRO_IDMAP:
    type: File
    inputBinding:
      position: 2

outputs:
  GAF_JSON:
    type: File
    outputBinding:
      glob: gaf.json
