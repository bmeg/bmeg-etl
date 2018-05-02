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


CURRENT_SHEET                    = 0
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

def CidToSid(name):
	return getresult("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/sids/TXT" % name)
def NameToCid(name):
	return getresult("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/TXT" % name)

def CidtoSMILES(name):
	return getresult("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/CanonicalSMILES/TXT" % name)





def clean(data):
	#print type(data)
	if type(data) == float:
		return str(int(data)).strip()
	else:
		return str(data).strip()

def getBioAssays(file):
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

	file = open("Bioassays.txt", "w")

	i = 0
	for cid in cids_array:
		cid = clean(cid)
		print "num: ", type(clean(num_array[i]))
		print "drug: ", type(drugs_array[i]), drugs_array[i]
		print "cid: ",type(cid)
		if i != 0:
			if cid.lower() != "error":
				if len(cid.split()) == 1:
					file.write(clean(num_array[i]) + ";" + cid + ";" + "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/assaysummary/CSV\n".format(cid))
				else:
					file.write(clean(num_array[i]) +  ";" + cid + ";ERROR\n")
			else:
				file.write(clean(num_array[i]) + ";"  ";" + cid + ";ERROR\n")




		i += 1
	file.close()

def Cid2Sids(file):
	"""https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/10000/sids/txt"""

	# Load in the workbook
	book = xlrd.open_workbook(file)

	#print book.sheet_names()
	sheet = book.sheet_by_index(CURRENT_SHEET)
	#print sheet.name
	drugs_array = sheet.col_values(COLUMN_CONTAINING_NAMES)
	cids_array = sheet.col_values(COLUMN_CONTAINING_CIDS)
	num_array = sheet.col_values(COLUMN_NUM)

	file = open("SIDs.txt", "w")

	i = 0
	for cid in cids_array:
		cid = clean(cid)
		print "num: ", clean(num_array[i])
		print "drug: ", drugs_array[i]
		if i != 0:
			if cid.lower() != "error":
				if len(cid.split()) == 1:
					file.write(clean(num_array[i]) +  ";" + cid + ";" + str(CidToSid(cid).split("\n"))+"\n")
				else:
					file.write(clean(num_array[i]) + ";" + cid + ";ERROR\n")
			else:
				file.write(clean(num_array[i]) + ";" + cid + ";ERROR\n")




		i += 1
	file.close()

def getCIDS(file):
	"""https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/10000/sids/txt"""

	# Load in the workbook
	book = xlrd.open_workbook(file)

	#print book.sheet_names()
	sheet = book.sheet_by_index(CURRENT_SHEET)
	#print sheet.name
	drugs_array = sheet.col_values(COLUMN_CONTAINING_NAMES)
	cids_array = sheet.col_values(COLUMN_CONTAINING_CIDS)
	num_array = sheet.col_values(COLUMN_NUM)

	file = open("CIDs.txt", "w")

	i = 0
	for drug in drugs_array:
		drug = clean(drug)
		print "num: ", clean(num_array[i])
		print "drug: ", drugs_array[i]
		if i != 0:
			if drug != None:
				try:
					val = NameToCid(drug)
					val2 = val.split("\n")
					print "cid: ", val2
					if len(val2) == 1:
						file.write(clean(num_array[i]) +  ";" + drug + ";" + str(val)+"\n")
					else:
						file.write(clean(num_array[i]) +  ";" + drug + ";" + str(val2)+"\n")

				except:
					file.write(clean(num_array[i]) + ";" + drug + ";Error\n")



			else:
				file.write(clean(num_array[i]) + ";" + cid + ";Error\n")




		i += 1
	file.close()


def getSMILES(file):

	# Load in the workbook
	book = xlrd.open_workbook(file)

	#print book.sheet_names()
	sheet = book.sheet_by_index(CURRENT_SHEET)
	#print sheet.name
	drugs_array = sheet.col_values(COLUMN_CONTAINING_NAMES)
	cids_array = sheet.col_values(COLUMN_CONTAINING_CIDS)
	num_array = sheet.col_values(COLUMN_NUM)

	file = open("SMILES.txt", "w")

	i = 0
	for cid in cids_array:
		cid = clean(cid)
		print "num: ", clean(num_array[i])
		print "drug: ", drugs_array[i]
		if i != 0:
			if cid.lower() != "error":
				if len(cid.split()) == 1:
					smiles = CidtoSMILES(cid)
					print smiles
					file.write(clean(num_array[i]) +  ";" + cid + ";" + str(smiles) + "\n")
				else:
					file.write(clean(num_array[i]) + ";" + cid + ";ERROR\n")
			else:
				file.write(clean(num_array[i]) + ";" + cid + ";ERROR\n")

		i += 1
	file.close()



def getBioRolesAndApplications(file):
	# Load in the workbook
	book = xlrd.open_workbook(file)

	#print book.sheet_names()
	sheet = book.sheet_by_index(CURRENT_SHEET)
	#print sheet.name
	drugs_array = sheet.col_values(COLUMN_CONTAINING_NAMES)
	cids_array = sheet.col_values(COLUMN_CONTAINING_CIDS)
	num_array = sheet.col_values(COLUMN_NUM)

	file = open("BioRoles.txt", "w")

	i = 0
	for cid in cids_array:
		cid = clean(cid)
		print "num: ", clean(num_array[i])
		print "drug: ", drugs_array[i]
		if i != 0:
			if cid.lower() != "error":
				if len(cid.split()) == 1:
					d = parseBioRolesAndApplications(cid)
					print d
					if len(d) == 5:
						file.write(clean(num_array[i]) +  ";" + cid + ";" + str(d[0]) + ";" + str(d[1]) + ";" + str(d[2]) + ";" + str(d[3]) + ";" + str(d[4]) + "\n")
					else:
						file.write(clean(num_array[i]) + ";" + cid + ";ERROR;ERROR;ERROR;ERROR;ERROR\n")
				else:
					file.write(clean(num_array[i]) + ";" + cid + ";ERROR;ERROR;ERROR;ERROR;ERROR\n")
			else:
				file.write(clean(num_array[i]) + ";" + cid + ";ERROR;ERROR;ERROR;ERROR;ERROR\n")

		i += 1
	file.close()



def parseBioRolesAndApplications(cid):
	

	ChEBI_ID = "ERROR"
	try:
		# Getting the ChEBI ID
		j = CIDData(cid)
		f = j.split("CHEBI:", 1)[1]
		ChEBI_ID = f.split("\"", 1)[0]
		ChEBI_website = getresult("http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:%s" % str(ChEBI_ID))
	except:
		d = ['ERROR','ERROR','ERROR','ERROR','ERROR']

	# Getting the info from ChEBI Oncology 
	try:
		biological_role = ChEBI_website.split("Biological Role", 1)[1].split("chebiId", 1)[1].split(">",1)[1].split("<",1)[0]
	except:
		biological_role = 'ERROR'

	try:
		biological_role_info = ChEBI_website.split("Biological Role", 1)[1].split("chebiId", 1)[1].split("\"roleDefinition\">",1)[1].split("<",1)[0]
	except:
		biological_role_info = 'ERROR'
	try:
		application = ChEBI_website.split("Biological Role", 1)[1].split("chebiId", 1)[1].split("chebiId=", 2)[2].split(">",1)[1].split("<",1)[0]
		if application.lower() == "application":
			try:
				application = ChEBI_website.split("Biological Role", 1)[1].split("chebiId", 1)[1].split("chebiId=", 2)[1].split(">",1)[1].split("<",1)[0]
			except:
				application = 'ERROR'
	except:
		application = 'ERROR'
	try:
		application_info = ChEBI_website.split("Biological Role", 1)[1].split("chebiId", 1)[1].split("chebiId=", 2)[2].split("\"roleDefinition\">",1)[1].split("<",1)[0]
	except:
		application_info = 'ERROR'
	d = [ChEBI_ID, biological_role, biological_role_info, application, application_info]

	####### Testing
	#print "ChEBI: ", ChEBI_ID
	#print "Biological Role: ", biological_role
	#print "Biological Role Info: ",biological_role_info
	#print "Application: ",application
	#print "Application Info: ",application_info


	return d



def CIDData(cid):
	return getresult("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON/" % str(cid))



if __name__ == "__main__":

	#Cid2Sids( "Drug_Conversion_test.xlsx")	
	#getBioAssays( "Drug_Conversion_MASTER.xlsx")
	#getCIDS("Drug_Conversion_test.xlsx")	
	#getSMILES("Drug_Conversion_MASTER.xlsx")
	getBioRolesAndApplications("../spreadsheets/Drug_Conversion_MASTER_UPDATED.xlsx")
	#parseBioRolesAndApplications("5790")
