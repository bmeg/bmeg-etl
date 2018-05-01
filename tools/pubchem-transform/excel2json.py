#!/usr/bin/python

import xlrd
from collections import OrderedDict
import simplejson as json
 
# Open the workbook and select the first worksheet
wb = xlrd.open_workbook('Drug_Conversion_MASTER.xlsx')
sh = wb.sheet_by_index(0)
 
# List to hold dictionaries
drugs_list = []

f = open("data.txt", "w")
 
# Iterate through each row in worksheet and fetch values into dict
for rownum in range(1, sh.nrows):
    drugs_list = []
    drugs = OrderedDict()
    row_values = sh.row_values(rownum)

    drugs['id']           = "compound:"+str(row_values[3])
    drugs['num']          = int(row_values[0])
    drugs['source']       = row_values[1]
    drugs['name']         = row_values[2]
    drugs['pubchem_cid']  = row_values[3]
    drugs['pubchem_sid']  = row_values[4]
    drugs['toxicity']     = row_values[5]
    drugs['bioassays']    = row_values[6]
    drugs['synonyms']     = "{}, {}, {}".format(row_values[7], row_values[8], row_values[9])
    drugs['broad_cpd_id'] = row_values[7]
    drugs['other_name']   = row_values[8]
    drugs['GNF_REG_ID']   = row_values[9]
    drugs['drug_id']      = row_values[10]
    drugs['smiles']       = row_values[11]
 
    #drugs_list.append(drugs)
 
    # Serialize the list of dicts to JSON
    j = json.dumps(drugs)
 
    # Write to file
    f.write(j + "\n")
f.close()