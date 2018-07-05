
##################
### CCLE

class: Workflow
cwlVersion: v1.0

$namespaces:
  bmeg: http://bmeg.io

inputs:
  DRUG_DATA:
    type: File
    bmeg:key: source/CCLE_NP24.2009_Drug_data_2015.02.24.csv
  PUBCHEM:
    type: File
    bmeg:key: source/ccle_pubchem.txt
  SAMPLE_TSV:
    type: File
    bmeg:key: source/CCLE_sample_info_file_2012-10-18.txt
  RPKM_GCT:
    type: File
    bmeg:key: source/CCLE_RNAseq_081117.rpkm.gct
  CCLE_MAF:
    type: File
    bmeg:key: source/ccle2maf_081117.txt

steps:
  ccle-drug-transform:
    run: ccle/ccle-drug-transform.cwl
    in:
      DRUG: DRUG_DATA
      #PUBCHEM: PUBCHEM
    out:
      - OUTPUT

  ccle-sample-transform:
    run: ccle/ccle-sample-transform.cwl
    in:
      SAMPLE_TSV: SAMPLE_TSV
    out:
      - OUTPUT

  ccle-expression-transform:
    run: gct/gct-transform.cwl
    in:
      SOURCE:
        default: ccle
      GCT: RPKM_GCT
      SCALE:
        default: RPKM
    out:
      - OUTPUT

  ccle-variant-transform:
    run: variant/maf-transform.cwl
    in:
      GENE_COLUMN:
        default: Hugo_Symbol
      SOURCE:
        default: ccle
      METHOD:
        default: CCLE
      BIOPREFIX:
        default: ccle
      CALLSETPREFIX:
        default: ccle
      MAF: CCLE_MAF
    out:
      - VARIANT
      - VARIANT_ANNOTATION
      - CALLSET

outputs:
  CCLE_VARIANT:
    type: File
    outputSource: ccle-variant-transform/VARIANT
