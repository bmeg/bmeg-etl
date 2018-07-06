
cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: bmeg/bmeg-etl:latest

baseCommand:
  - python
  - /opt/ga4gh-variant.py

arguments:
  - "--multi"
  - maf_data

inputs:
  MAF:
    type: File
    inputBinding:
      prefix: "--maf"

  BIOPREFIX:
    type: string
    inputBinding:
      prefix: "--bioPrefix"

  CALLSETPREFIX:
    type: string
    inputBinding:
      prefix: "--callSetPrefix"

  GENECOL:
    type: string?
    inputBinding:
      prefix: "--gene"

  CENTER_COLUMN:
    type: string?
    inputBinding:
      prefix: "--center"

  SOURCE:
    type: string?
    inputBinding:
      prefix: "--source"

  GENOME_BUILD:
    type: string?
    inputBinding:
      prefix: "--genome"

  METHOD:
    type: string?
    inputBinding:
      prefix: "--method"
  GZIP:
    type: boolean?
    inputBinding:
      prefix: "--gz"

outputs:
  VARIANT:
    type: File
    outputBinding:
      glob: "maf_data.bmeg.Variant.json"
  VARIANT_ANNOTATION:
    type: File
    outputBinding:
      glob: "maf_data.bmeg.VariantAnnotation.json"
  CALLSET:
    type: File
    outputBinding:
       glob: "maf_data.bmeg.CallSet.json"
