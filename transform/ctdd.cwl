


##################
### CCLE

class: Workflow
cwlVersion: v1.0

$namespaces:
  bmeg: http://bmeg.io

inputs:
  CTRP_ZIP:
    type: File
    bmeg:key: source/CTRPv2.0_2015_ctd2_ExpandedDataset.zip

steps:

  ctdd-response-transform:
    run: ctdd/ctdd-transform.cwl
    in:
      CTRP_ZIP: CTRP_ZIP
    out:
      - RESPONSE

outputs:
  RESPONSE_JSON:
    type: File
    outputSource: ctdd-response-transform/RESPONSE
