#!/bin/bash

# pip install protobuf
# pip install pandas
# pip install xlrd

# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Cell_Lines_Details.xlsx
# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/v17_fitted_dose_response.xlsx
# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/GDSC-CCLE-CTRP_conversion.xlsx
# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Screened_Compounds.xlsx
# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/v17a_public_raw_data.xlsx

./gdsc-convert.py GDSC-CCLE-CTRP_conversion.xlsx Cell_Lines_Details.xlsx Screened_Compounds.xlsx v17a_public_raw_data.xlsx v17_fitted_dose_response.xlsx gdsc_pubchem.table
