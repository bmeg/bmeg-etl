

cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: gdsc-transform:latest

baseCommand:
  - /opt/gdsc-convert.py

arguments:
  - position: 6
    valueFrom: /opt/gdsc_pubchem.table

inputs:
  CONVERSION: # GDSC-CCLE-CTRP_conversion.xlsx
    type: File
    inputBinding: 
      position: 1

  CELL_LINES: # Cell_Lines_Details.xlsx
    type: File
    inputBinding: 
      position: 2
  
  COMPOUNDS: # Screened_Compounds.xlsx
    type: File
    inputBinding: 
      position: 3
  
  RAW_DATA: # v17a_public_raw_data.xlsx 
    type: File
    inputBinding:
      position: 4
      
  RESPONSE: # v17_fitted_dose_response.xlsx
    type: File
    inputBinding:
      position: 5
    
  
outputs:
  GDSC_JSON:
    type: 
      type: array
      items: File
    outputBinding:
      glob: "*.json"
