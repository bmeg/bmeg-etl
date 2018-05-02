#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/gdc-transform:latest

baseCommand:
  - python
  - /opt/convert-clinical.py

arguments:
  - "--out"
  - ./
  - "--pubchem"
  - "/opt/tcga_pubchem.map"

inputs:
  TARBALL:
    type: File
    inputBinding:
      prefix: "--tar"

outputs:
  INDIVIDUAL:
    type: File
    outputBinding:
      glob: tcga.Individual.json

  CLINICAL_FOLLOWUP:
    type: File
    outputBinding:
      glob: tcga.ClinicalFollowup.json

  DRUG_THERAPY:
    type: File
    outputBinding:
      glob: tcga.DrugTherapy.json

  RADIATION_THERAPY:
    type: File
    outputBinding:
      glob: tcga.RadiationTherapy.json
