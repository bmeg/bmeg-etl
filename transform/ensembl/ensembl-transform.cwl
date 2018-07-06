
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest


baseCommand:
  - go
  - run
  - /opt/run.go

inputs:
  GAF_GZ:
    type: File
    inputBinding:
      position: 1

outputs:
  TRANSCRIPT:
    type: File
    outputBinding:
      glob: Transcript.json
  GENE:
    type: File
    outputBinding:
      glob: Gene.json
  EXON:
    type: File
    outputBinding:
      glob: Exon.json
