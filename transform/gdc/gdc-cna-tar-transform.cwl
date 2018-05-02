#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: biostream/gdc-transform:latest


baseCommand: 
  - python
  - /opt/convert-seg-tar.py

inputs:  
  FILE_MAP:
    type: File
    inputBinding:
      prefix: "--filemap"
  
  TARBALL:
    type: File
    inputBinding:
      prefix: "--tar"
  
  GTF:
    type: File
    inputBinding:
      prefix: "--gtf"  

outputs:
  CNASegment:
    type: File
    outputBinding:
      glob: out.bmeg.CNASegment.json
  CNACallSet:
    type: File
    outputBinding:
      glob: out.bmeg.CNACallSet.json