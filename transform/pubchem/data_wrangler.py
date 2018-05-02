#!/usr/bin/python2

"""
data_wrangler.py

Author:
    Theodore J. LaGrow

Usage:
    Attempt to extract data and metadata from PubChem and
    various other sources.
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
import threading
import pandas as pd
import time


# GLOBAL VARIABLES
MULTITHREAD = True
EXCELPATH = "../spreadsheets/Drug_Conversion_MASTER_UPDATED.xlsx"
SHEETNAME = "CTDD"
THREADGROUPNUM = 5

#############################################################
# Check to see if user wants single or multithreaded

def start():
    """ init """
    global MULTITHREAD, THREADGROUPNUM
    ans1 = raw_input("Do you want to multithread [yes/no]? ")
    if ans1.lower() == "yes" or ans1.lower() == "y":
        ans2 = raw_input("How many threads do you want at a time? (recommended <5) ")
        try: 
            assert int(ans2)
            THREADGROUPNUM = int(ans2.strip())            
        except:
            print("Please input an integer. Try again.")
            start()

    elif ans1.lower() == "no" or ans1.lower() == "n":
        MULTITHREAD = False
    else:
        print("Please use either 'yes' or 'no'")
        start()


#############################################################
# Formating functions 

def getresult(url):
    """ Simple function to attempt to connect to a website """
    try:
        connection = urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        return "Error"
    else:
        return connection.read().rstrip()

def OpenExcelWorkbookSheet(path, name_of_sheet):
    """ Returns the dataframe of the excel sheet """
    return pd.DataFrame(pd.read_excel(path, name_of_sheet))

def pausing(n):
    if n % THREADGROUPNUM == 0:{time.sleep(1)}

def checkMulti(output):
    if len(output.split("\n")) != 1:
        return output.split("\n")
    else:
        return output.rstrip()

def StripNames(name):
    """ 
    input: drug name for an input
    output: a name stripped of anything in parenthesis, these are typically notes from the person
        who inputed the data 
    use: Mostly for G2P 
    """
    name = str(name)
    name = re.sub(r'\([^)]*\)', '', name)
    #name = re.sub('inhibitor', '', name) #specific word someone put into the data that messes with the input for api and scrapping
    #name = re.sub('inhibitors', '', name) #specific word someone put into the data that messes with the input for api and scrapping
    return name


#############################################################
# Information retreiving functions

def NameToCid(name):
    """  
    input: name of drug
    output: CID from PubChem using API
    """
    return getresult("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/TXT" % name)

def NameToSid(name):
    """
    input: name of drug
    output: SID from PubChem using API
    """
    return getresult("http://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/%s/sids/TXT" % name)


#############################################################


def getCID(id, name, file):
    cid = checkMulti(NameToCid(StripNames(name)))
    file.write("{}\t{}\t{}\n".format(id, name, cid))
    print "Thread %s %s\n" % (id, cid)
    return

def getSID(id, name, file):
    sid = checkMulti(NameToSid(StripNames(name)))
    file.write("{}\t{}\t{}\n".format(id, name, sid))
    print "Thread %s %s\n" % (id, sid)
    return



if __name__ == "__main__":

    start()

    ##############################################################
    # Please take a minute and check to make sure you are running
    # the functions you want
    ##############################################################

    sheet = OpenExcelWorkbookSheet(EXCELPATH, SHEETNAME)

    file =  open("Output.tsv", "w")
    i = 1
    for n in sheet.Name:
        if MULTITHREAD == True:
            pausing(i)
            cid_thread = threading.Thread(target=getSID, args=(i,n,file,))
            cid_thread.start()
            i+=1
        else:
            getCID(i,n,file)
            i+=1

    time.sleep(4)
    file.close()

