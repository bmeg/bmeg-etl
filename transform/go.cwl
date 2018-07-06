##################
### Gene Ontology

class: Workflow
cwlVersion: v1.0

$namespaces:
  bmeg: http://bmeg.io

inputs:
  OBO:
    type: File
    bmeg:key: source/go.obo
  GAF:
    type: File
    bmeg:key: source/goa_human.gaf.gz
  UNIPRO_MAP:
    type: File
    bmeg:key: source/HUMAN_9606_idmapping.dat.gz

steps:
  go-transform:
    run: go/go-transform.cwl
    in:
      OBO: OBO
    out:
      - GO_JSON

  gaf-transform:
    run: go/gaf-transform.cwl
    in:
      GAF: GAF
      UNIPRO_IDMAP: UNIPRO_MAP
    out:
      - GAF_JSON

outputs:
  GOPROTO:
    type: File
    outputSource: go-transform/GO_JSON
  GAF_JSON:
    type: File
    outputSource: gaf-transform/GAF_JSON
