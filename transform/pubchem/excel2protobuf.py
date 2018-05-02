#!/usr/bin/python

import xlrd
import simplejson as json
import phenotype_pb2
import sys
from google.protobuf.json_format import MessageToJson




def CheckUnicode(s):
	if isinstance(s, float):
		s = str(int(s))
	if isinstance(s, unicode):
		return s
	else:
		return unicode(s, "UTF-8")

# Main procedure:  Reads the data collected from a file,
#   construct a message row by row, then write it to a text file.

if len(sys.argv) != 2:
  print "Usage:", sys.argv[0], "BINARY_DATA_FILE"
  sys.exit(-1)

compound_data = phenotype_pb2.Compound()

# Read the existing file
try:
  f = open(sys.argv[1], "rb")
  compound_data.ParseFromString(f.read())
  f.close()
except IOError:
  print sys.argv[1] + ": Could not open file.  Creating a new one."


# Open the workbook and select the first worksheet
wb = xlrd.open_workbook('../spreadsheets/Drug_Conversion_MASTER_UPDATED.xlsx')
sh = wb.sheet_by_index(0)

# Iterate through each row in worksheet and fetch values into dict
for rownum in range(1, sh.nrows):

    compound_data = phenotype_pb2.Compound() #construct new message
    row_values = sh.row_values(rownum)
    
    compound_data.gid         = CheckUnicode("compound:"+str(row_values[3]))
    compound_data.id          = CheckUnicode(str(int(row_values[0])))
    compound_data.source      = CheckUnicode(row_values[1])
    compound_data.name        = CheckUnicode(row_values[2])
    compound_data.pubchem_cid = CheckUnicode(row_values[3])
    compound_data.pubchem_sid = CheckUnicode(row_values[4])
    compound_data.toxicity    = CheckUnicode(row_values[6])
    compound_data.bioassays   = CheckUnicode(row_values[7])
    compound_data.smiles      = CheckUnicode(row_values[12])

    #adding fields to synonyms: broad_cpd_id = 8, other_name = 9, GNF_REG_ID = 10, drug_id = 11
    for i in range(8, 12):
    	if row_values[i] != "":
    		compound_data.synonyms.append(CheckUnicode(row_values[i]))


    # Write the messages to binary
    f1 = open(sys.argv[1], "a")
    f1.write(compound_data.SerializeToString())
    f1.close()

    # Generate the json format of the messages
    jsonObj = json.dumps(json.loads(MessageToJson(compound_data)))


    ####### TESTING
    #if rownum == 462:
    #	print "*******************yes"
    #	print row_values[8]
    #	print type(row_values[8])
    #	print jsonObj

    # Wrine the messages in JSON by line
    f2 = open("JSONprotobuf.txt", "a")
    f2.write(jsonObj + "\n")
    f2.close()


    print rownum #printing for sanity