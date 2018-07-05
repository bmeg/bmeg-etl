cwlVersion: v1.0
class: Workflow

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  PROJECT_ID:
    type: string
  GTF:
    type: File

steps:
  download-gdc-expression:
    in:
      DATA_TYPE:
        default: fpkm
      PROJECT_ID: PROJECT_ID
      OUTPUT_NAME:
        valueFrom: $(inputs.PROJECT_ID + ".fpkm.tar.gz")
    out:
      - ARCHIVE
      - FILEMAP
    run: src/github.com/biostream/gdc-extract/gdc-bulk-download.cwl

  download-gdc-biosample:
    in:
      DATA_TYPE:
        default: biospecimen
      PROJECT_ID: PROJECT_ID
      OUTPUT_NAME:
        valueFrom: $(inputs.PROJECT_ID + ".biosample.tar.gz")
    out:
      - ARCHIVE
      - FILEMAP
    run: src/github.com/biostream/gdc-extract/gdc-bulk-download.cwl

  download-gdc-clinical:
    in:
      DATA_TYPE:
        default: clinical
      PROJECT_ID: PROJECT_ID
      OUTPUT_NAME:
        valueFrom: $(inputs.PROJECT_ID + ".clinical.tar.gz")
    out:
      - ARCHIVE
      - FILEMAP
    run: src/github.com/biostream/gdc-extract/gdc-bulk-download.cwl

  gdc-rna-transform:
    run: src/github.com/biostream/gdc-transform/gdc-expression-tar-transform.cwl
    in:
      TARBALL: download-gdc-expression/ARCHIVE
      FILE_MAP: download-gdc-expression/FILEMAP
    out:
      - DATA

  gdc-clinical-transform:
    run: src/github.com/biostream/gdc-transform/gdc-clinical-tar-transform.cwl
    in:
      TARBALL: download-gdc-clinical/ARCHIVE
    out:
      - INDIVIDUAL
      - CLINICAL_FOLLOWUP
      - DRUG_THERAPY
      - RADIATION_THERAPY



outputs:
  EXPRESSION_TAR:
    type: File
    outputSource: download-gdc-expression/ARCHIVE
  INDIVIDUAL:
    type: File
    outputSource: gdc-clinical-transform/INDIVIDUAL
  CLINICAL_FOLLOWUP:
    type: File
    outputSource: gdc-clinical-transform/CLINICAL_FOLLOWUP
  DRUG_THERAPY:
    type: File
    outputSource: gdc-clinical-transform/DRUG_THERAPY
  RADIATION_THERAPY:
    type: File
    outputSource: gdc-clinical-transform/RADIATION_THERAPY
  FPKMS:
    type: File
    outputSource: gdc-rna-transform/DATA
