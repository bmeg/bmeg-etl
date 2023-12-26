
cwlVersion: v1.2

requirements:
  - class: DockerRequirement
    dockerImageId: pharmaco-gx:4.3.0
  - class: NetworkAccess
    networkAccess: true

class: CommandLineTool
baseCommand: [Rscript, /usr/src/download_pharmaco.R, "-s", "./"]

inputs: {}

outputs:
  profiles:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*"