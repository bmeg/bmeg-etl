
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest

baseCommand:
  - "python3.7"
  - "/opt/ensembl-transform.py"

outputs:
  GENE:
    type: File
    outputBinding:
      glob: ensembl.Gene.Vertex.json
  TRANSCRIPT:
    type: File
    outputBinding:
      glob: ensembl.Transcript.Vertex.json
  TRANSCRIPTFOR:
    type: File
    outputBinding:
      glob: ensembl.TranscriptFor.Edge.json
  EXON:
    type: File
    outputBinding:
      glob: ensembl.Exon.Vertex.json
  EXONFOR:
    type: File
    outputBinding:
      glob: ensembl.ExonFor.Edge.json
