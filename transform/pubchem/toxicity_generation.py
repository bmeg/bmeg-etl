#!/usr/bin/python

"""
toxicity_generation.py

Author:
Theodore J. LaGrow
lagrow@ohsu.edu

Usage:
Attempt to extract toxicity reports from compound ID 
and substance ID 

"""



import os
import sys
import urllib2
import xlrd
import re
from lxml import html
import requests
from bs4 import BeautifulSoup
from pprint import pprint
import json



# CTDD -> Sheet 1, Column for IDs-> 3

CURRENT_SHEET                    = 1
COLUMN_CONTAINING_NAMES          = 2
COLUMN_CONTAINING_CIDS           = 3
COLUMN_NUM                       = 0




def getresult(url):
	try:
		connection = urllib2.urlopen(url)
	except urllib2.HTTPError, e:
		return "Error"
	else:
		return connection.read().rstrip()

def NameToCid(name):
	return getresult("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/TXT" % name)

def NameToSid(name):
	return getresult("http://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/%s/sids/TXT" % name)

def GetToxicityDatafromPubChem(cid):
	return getresult("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON/" % cid)

def GetAssayDatafromPubChem(cid):
	#return getresult("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON" % cid)
	return getresult("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/CSV" % cid)

def tryBMEG(gid):
	return getresult("https://bmeg.compbio.ohsu.edu/?gid=gene:%s" % gid)


def clean(data):
	#print type(data)
	if type(data) == float:
		return str(int(data)).strip()
	else:
		return str(data).strip()

def ToxicityReport(j):

	#j = GetToxicityDatafromPubChem(j)
	data = json.loads(j)
	field = data[u'Record'][u'Section']
	for n in range(len(field)):
		#print pprint(field[n])

		#print field[n][u'Description']
		f = field[n][u'Description'].split()
		if f[0].strip() == "Toxicity":
			return field[n][u'Section']



def BioAssaysReport(j):
	""" Not quite, the data says to go external """
	""" found it... """
	""" pretty much useless now... """
	j = GetToxicityDatafromPubChem(j)
	data = json.loads(j)
	field = data[u'Record'][u'Section']
	for n in range(len(field)):
		f = field[n][u'Description'].split()
		for word in f:
			if word.lower() == "bioassay":
				return field[n]



def printLineBreaks():
	print "\n\n\n\n\n\n**************************\n\n\n\n\n\n"

def parse(file):
	og_errors = 0
	new_errors = 0
	actual_tox_reports = 0

	# Load in the workbook
	book = xlrd.open_workbook(file)

	#print book.sheet_names()
	sheet = book.sheet_by_index(CURRENT_SHEET)
	#print sheet.name
	drugs_array = sheet.col_values(COLUMN_CONTAINING_NAMES)
	cids_array = sheet.col_values(COLUMN_CONTAINING_CIDS)
	num_array = sheet.col_values(COLUMN_NUM)
	file = open("toxicity.txt", "w")


	i = 0
	for cid in cids_array:
		cid = clean(cid)
		if i != 0:
			#print "Drug testing: ", drugs_array[i]
			#print "CID tested: ", cid.lower()
			if str(cid).lower() != "error":
				try:
					j = GetToxicityDatafromPubChem(cid) #to get the JSON data once
					#print ToxicityReport(j)
					file.write(clean(num_array[i]) + "@"  + cid + "@" + str(ToxicityReport(j))+"\n")


					#printLineBreaks()
					actual_tox_reports += 1
				except ValueError:
					#print "Error in accessing Toxicity"
					file.write(clean(num_array[i]) + "@"  + cid + "@Error in accessing Toxicity\n")

					#printLineBreaks()
					new_errors += 1
			else:
				#print "Error in Origional Data"
				file.write(clean(num_array[i]) + "@"  + cid + "@Error in Origional Data\n")

				#printLineBreaks()
				og_errors += 1
		print i, "\n" 
		i += 1

	print
	print
	print "OG  Errors: ", og_errors
	print "New Errors: ", new_errors
	print "Act Tox Re: ", actual_tox_reports
	file.close()

	return 0




if __name__ == "__main__":

	#print pprint(ToxicityReport("31703"))
	#print pprint( BioAssaysReport("7326481") )
	#print GetAssayDatafromPubChem("7326481")
	#parse( "Drug_Conversion_MASTER.xlsx")


