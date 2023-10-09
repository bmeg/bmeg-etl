
cwlVersion: v1.2
class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerImageId: pharmaco-gx:4.3.0

baseCommand: [Rscript, /usr/src/drug_response_single.R,  --out=./]

inputs:
  input:
    type: File
    inputBinding:
      prefix: --input=
      separate: false

outputs:
  response:
      type: File
      outputBinding:
        glob: response.tsv
  samples:
      type: File
      outputBinding:
        glob: samples.tsv
  treatments:
      type: File
      outputBinding:
        glob: treatments.tsv
