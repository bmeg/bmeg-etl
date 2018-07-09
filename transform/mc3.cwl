##################
### MC3

class: Workflow
cwlVersion: v1.0

$namespaces:
  bmeg: http://bmeg.io

inputs:
  MC3_MAF:
    type: File
    bmeg:key: source/mc3.v0.2.8.PUBLIC.maf.gz

steps:

  mc3-transform:
    run: variant/maf-transform.cwl
    in:
      CENTER_COLUMN:
        default: "CENTERS"
      GENE_COLUMN:
        default: "Gene"
      SOURCE:
        default: tcga
      GZIP:
        default: true
      BIOPREFIX:
        default: "tcga"
      CALLSETPREFIX:
        default: mc3
      MAF: MC3_MAF
    out:
      - VARIANT
      - VARIANT_ANNOTATION
      - CALLSET

outputs:
  VARIANT_JSON:
    type: File
    outputSource: mc3-transform/VARIANT
    bmeg:key: biostream/mc3/mc3.Variant.json
  VARIANT_ANNOTATION_JSON:
    type: File
    outputSource: mc3-transform/VARIANT_ANNOTATION
    bmeg:key: biostream/mc3/mc3.VariantAnnotation.json
  CALLSET_JSON:
    type: File
    outputSource: mc3-transform/CALLSET
    bmeg:key: biostream/mc3/mc3.CallSet.json
