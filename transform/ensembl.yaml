##################
### ENSEMBL

class: Workflow
cwlVersion: v1.0


$namespaces:
  bmeg: http://bmeg.io

inputs: []

steps:
  ensembl-extract:
    run: curl/curl.cwl
    in:
      URL:
        default: ftp://ftp.ensembl.org/pub/release-90/gff3/homo_sapiens/Homo_sapiens.GRCh38.90.chr.gff3.gz
      NAME:
        default: Homo_sapiens.GRCh38.90.chr.gff3.gz
    out:
      - OUTPUT

  ensembl-transform:
    run: ensembl/ensembl-transform.cwl
    in:
      GAF_GZ: ensembl-extract/OUTPUT
    out:
      - TRANSCRIPT
      - GENE
      - EXON

outputs:
  TRANSCRIPT:
    type: File
    outputSource: ensembl-transform/TRANSCRIPT
  GENE:
    type: File
    outputSource: ensembl-transform/GENE
  EXON:
    type: File
    outputSource: ensembl-transform/EXON
