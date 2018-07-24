
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest
  
baseCommand:
  - /opt/gdc-scan-cases.py

stdout: gdc-cases.json

outputs:
  CASES:
    type: stdout
