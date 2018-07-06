
##################
### GTEX

class: Workflow
cwlVersion: v1.0

$namespaces:
  bmeg: http://bmeg.io

inputs:
  BIOBANK:
    type: File
    bmeg:key: source/gtex/biobank_collection_20180116_031101.txt
  GTEX:
    type: File
    bmeg:key: source/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz

steps:
  gtex-sample-transform:
    run: gtex/gtex-transform.cwl
    in:
      BIOBANK: BIOBANK
    out:
      - INDIVIDUAL
      - BIOSAMPLE

  gtex-expression-transform:
    run: gct/gct-transform.cwl
    in:
      SOURCE:
        default: gtex
      GZIPPED:
        default: true
      GCT: GTEX
      SCALE:
        default: RPKM
    out:
      - OUTPUT

outputs:
  GeneExpression:
    type: File
    outputSource: gtex-expression-transform/OUTPUT
    bmeg:key: biostream/gtex/gtex.GeneExpression.json
  Biosample:
    type: File
    outputSource: gtex-sample-transform/BIOSAMPLE
    bmeg:key: biostream/gtex/gtex.Biosample.json
  Individual:
    type: File
    outputSource: gtex-sample-transform/INDIVIDUAL
    bmeg:key: biostream/gtex/gtex.Individual.json
