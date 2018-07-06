#!/usr/bin/env python

import json, sys, argparse, os
import csv #for drug data
import string
import re
import pandas
####### Added by TL, 7/7
import xlrd #extracting data from the excel spreadsheet
import math #used checking for NaN
#######
import numpy as np

def process(path, sheet):
    print "Reading in the data is stupid, thanks Obama...\n\n"


    data_xls = pandas.read_excel(path, sheet, index_col=None)

    data_xls.to_csv('csvfile1.csv', encoding='utf-8', index=False, sep="\t")

    exceldata_df = pandas.read_table('csvfile1.csv')
    print "\n\nDataframe ready"

    file = open("outputMaster.tsv", "w")



    #i = 0
    old_row = ""
    col = exceldata_df.columns.values
    print "Starting processing"

    for index, row in exceldata_df.iterrows():
    	if index == 1804:
    		exit("Finished!")
    	if index != 0:
    		if row.Repeat_Prev == "Same":
    			for j in range(len(row)):
    				if row[j] != old_row[j]:
    					#print "Row: ", row[j]
    					#print "Old: ", old_row[j]
    					try:
    						assert np.isnan(row[j])              #row == 1
    						try: 
    							assert np.isnan(old_row[j])      #old_row == 1
    							row[j] = unicode('', "UTF-8")
    						except:                              #old_row == 0
    							row[j] = old_row[j]
    					except:                                  #row == 0
    						try: 
    							assert np.isnan(old_row[j])      #old_row == 1
    							row[j] = row[j]
    						except:                              #old_row == 0
    							row[j] = [row[j], old_row[j]]
    				else:
    					row[j]
    				#print "ROW[J]: ", row[j]
    		r=""
    		for item in row:
    			r = r + str(item) + "\t"
    		file.write(str(r))
    		file.write("\n")
    		old_row = row
    			


    	else:
    		old_row = row #starting main loop
    	#i += 1


    file.close()






if __name__ == "__main__":

	process("../spreadsheets/Drug_Conversion_MASTER_UPDATED.xlsx", "Master_noerror")
