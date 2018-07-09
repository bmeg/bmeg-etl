##################
### ENSEMBL

class: Workflow
cwlVersion: v1.0


$namespaces:
  bmeg: http://bmeg.io

inputs: []

steps:
  pubmed-scan:
    run: pubmed/pubmed-list.cwl
    in: {}
    out:
      - DATA

outputs:
  PUBMED:
    type: File
    outputSource: pubmed-scan/DATA
    bmeg:key: source/pubmed.list
