

# http://www.gtexportal.org/static/datasets/biobank/downloads/biobank_collection_20170329_093753.txt
# http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz


cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/gtex-transform:latest

baseCommand:
  - /opt/gtex_convert.py
  - "--format"
  - json


arguments:
  - "--multi"
  - "$(inputs.OUTPATH)"

inputs:

  SAMPLES:
    type: File
    inputBinding:
      prefix: "--sample"
  OUTPATH:
    type: [string, "null"]
    default: "ccle-data"

outputs:
  OUTPUT:
    type: File[]
    outputBinding:
      glob: "*.json"