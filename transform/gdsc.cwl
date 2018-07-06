
class: Workflow
cwlVersion: v1.0

$namespaces:
  bmeg: http://bmeg.io

inputs:
  CONVERSION:
    type: File
    bmeg:key: source/gdsc/GDSC-CCLE-CTRP_conversion.xlsx
  CELL_LINES:
    type: File
    bmeg:key: source/gdsc/Cell_Lines_Details.xlsx
  COMPOUNDS:
    type: File
    bmeg:key: source/gdsc/Screened_Compounds.xlsx
  RAW_DATA:
    type: File
    bmeg:key: source/gdsc/v17.3_public_raw_data.csv
  RESPONSE:
    type: File
    bmeg:key: source/gdsc/v17.3_fitted_dose_response.xlsx

steps:
  gdsc-transform:
    run: gdsc/gdsc-transform.cwl
    in:
      CONVERSION: CONVERSION
      CELL_LINES: CELL_LINES
      COMPOUNDS: COMPOUNDS
      RAW_DATA: RAW_DATA
      RESPONSE: RESPONSE
    out:
      - GDSC_JSON

outputs:
  RESPONSE_JSON:
    type: File[]
    outputSource: gdsc-transform/GDSC_JSON
