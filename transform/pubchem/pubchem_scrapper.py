#!/usr/bin/python

"""
pubchem_scrapper.py

Author:
    Theodore J. LaGrow

Usage:
    Attempt to extract compound ID and substance ID 
    from name in given spreadsheet.  Completion of 
    data tool.

"""

import os
import xlrd
#from lxml import html
from bs4 import BeautifulSoup

# Global Variables 
#
# All of the sheets have different columns for the drug name
# so all of the sheet need to be changed according to the column name
# Please look at the data that you need to parse before running
#
# [sheet name]       ,    [column containing name]
# SMMART         (0) ,    0
# CTDD           (1) ,    0
# CCLE           (2) ,    1 (sometimes 2)
# SangerCancerRX (3) ,    1
# G2P            (4) ,    0

NAME_OF_FILE              = "Drug_Conversion_test.xlsx"
CURRENT_SHEET             = 1
COLUMN_CONTAINING_NAMES   = 1
DRUGS_THAT_OUTPUT_LISTS   = 0


def CheckBeforeRunning(file):
    """
    This function is to help with verifying what you want to parse before parsing to not 
    waste time, 

    pretty dumb function honestly, took waaaay too long create vs. the time it 
    takes to run.  One of the problems is that the data wrangling is comstantly shifting 
    and evolving so even the data file is different and doesn't use these hard coded variables
    """
    global NAME_OF_FILE, CURRENT_SHEET, COLUMN_CONTAINING_NAMES

    # Load in the workbook
    book = xlrd.open_workbook(file)
    sheet = book.sheet_by_index(CURRENT_SHEET)
    check = True

    while check == True:
        print "\nThe file being used is: ", file
        print "The sheet name is:      ", sheet.name, "\n"

        answer = raw_input("Is the file and sheet name correct? (yes/no) ")
        if answer.lower() == "yes":
            break
        elif answer.lower() == "no":
            while True:
                answer2 = raw_input("\nDo you know the path for the file? (yes/no) ")
                if answer2.lower() == "yes":
                    path = raw_input("Path: ")
                    NAME_OF_FILE = path
                    if os.path.exists(NAME_OF_FILE) == True:
                        print "\n*File path exists*\n"
                        while True:
                            array = []
                            file_sheet = xlrd.open_workbook(NAME_OF_FILE)
                            for n in range(file_sheet.nsheets):
                                array.append(file_sheet.sheet_by_index(n).name)

                            print "Number of sheets: ", file_sheet.nsheets, " ", array
                            print "\nPlease provide the sheet number"
                            sheet_number = raw_input("(starting at zero): ")
                            val = int(sheet_number)
                            
                            if val < file_sheet.nsheets:
                                confirmation = raw_input("\nIs %s the correct sheet? (yes/no) " % \
                                    file_sheet.sheet_by_index(val).name)
                                if confirmation.lower() == "yes":
                                    CURRENT_SHEET = val
                                    col_num = raw_input("\nWhat is the column containing drug names? ")
                                    COLUMN_CONTAINING_NAMES = int(col_num)
                                    check = False
                                    break

                                elif confirmation.lower() == "no":
                                    print "\nPlease try inputing the sheet number again\n"
                            else:
                                print "Error: inputed sheet number out of range, try again \n"

                    else:
                        print "Error: File does not exist. Please provide correct path"
                        NAME_OF_FILE = file
                    break

                elif answer2.lower() == "no":
                    exit("\nFind the path for your file and restart this program when you are ready\n")
                else:
                    print "please use 'yes' or 'no'"

        else:
            print "please use 'yes' or 'no'"
    print "\n***************************\nStarting data collection...\n***************************\n"


def ListFileSheetColumn():
    """ Use for testing """
    global NAME_OF_FILE, CURRENT_SHEET, COLUMN_CONTAINING_NAMES
    print "File: ", NAME_OF_FILE
    print "Sheet: ", xlrd.open_workbook(NAME_OF_FILE).sheet_by_index(CURRENT_SHEET).name
    print "Column: ", COLUMN_CONTAINING_NAMES
    print

def getresult(url):
    """ Simple function to attempt to connect to a website """
    try:
        connection = urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        return "Error"
    else:
        return connection.read().rstrip()

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

def Check4List(output):
    """
    input: drug contigent of having multiple drug in origional entry
    output: a list of the drugs from the origional entry, even if one drug from the origional entry
    """
    global DRUGS_THAT_OUTPUT_LISTS
    if len(output.split("\n")) != 1:
        DRUGS_THAT_OUTPUT_LISTS += 1
        return output.split("\n")
    else:
        return output





def ParseExcel(file):
    """
    This function returns a file named [sheetname]_ids_raw.txt 
    that has the name of the drug, if it is a compound or substance
    and if it returns the ID number.  The values are seperated 
    by commas to open in an spreadsheet to copy over easy.
    """

    # Load in the workbook
    book = xlrd.open_workbook(file)

    #print book.sheet_names() 

    sheet = book.sheet_by_index(CURRENT_SHEET)
    #print sheet.name

    drugs_array = sheet.col_values(COLUMN_CONTAINING_NAMES)
    num_of_drugs = len(drugs_array)

    file_name = sheet.name + "_ids_raw.txt"

    file = open(file_name, "w")

    error_counter = 0.0

    i = 0
    for drug in drugs_array:
    	og_drug = drug
        drug = StripNames(drug)
        drug_id_list = None
        drug_type_delineation_list = None

        if i != 0:
            drug_split = drug.split(" + ")
            #print "drug_split: ", drug_split
            num_of_drugs += (len(drug_split) - 1) #error checking, subtract the header cell in file
            for d in range(len(drug_split)):


                input_drug = drug_split[d].strip()

                cid_num = Check4List( NameToCid(input_drug) )
                cid_num = str(cid_num)
                sid_num = Check4List( NameToSid(input_drug) )
                sid_num = str(sid_num)

                
                if cid_num != "Error" and sid_num != "Error":
                    if drug_type_delineation_list == None:
                        drug_type_delineation_list = "CS"
                    else:
                        drug_type_delineation_list += ",CS"

                    if drug_id_list == None:
                        drug_id_list = cid_num + ";" + sid_num
                    else:
                        drug_id_list += "," + cid_num + ";" + sid_num 
                   
                     


                elif cid_num != "Error":
                    if drug_type_delineation_list == None:
                        drug_type_delineation_list = "C"
                    else:
                        drug_type_delineation_list += ",C"

                    if drug_id_list == None:
                        drug_id_list = cid_num
                    else:
                        drug_id_list += "," + cid_num

                elif sid_num != "Error":

                    if drug_type_delineation_list == None:
                        drug_type_delineation_list = "S"
                    else:
                        drug_type_delineation_list += ",S"

                    if drug_id_list == None:
                        drug_id_list = sid_num
                    else:
                        drug_id_list += "," + sid_num

                else:
                    error_counter += 1
                    if drug_type_delineation_list == None:
                        drug_type_delineation_list = "E"
                    else:
                        drug_type_delineation_list += ",E"

                    if drug_id_list == None:
                        drug_id_list = "Error"
                    else:
                        drug_id_list += ",Error"

            file.write(drug + ";" + drug_type_delineation_list + ";" + str(drug_id_list) + "\n")

            print og_drug, ":  ", drug_type_delineation_list, "  ", drug_id_list

        i += 1

    file.close()

    

    reliability = error_counter / num_of_drugs
    print "Errors occured: ", error_counter
    print "Number of drugs: ", num_of_drugs
    print "Percentage of errors: ", reliability
    print "Items to still go through... ", DRUGS_THAT_OUTPUT_LISTS
    

    # Uncomment to help generate report on errors occurring
    
    report_file = open("error_report.txt", "a")
    report_file.write( file_name                                                               +  "\n"     )
    report_file.write( "    Number of drugs: "      + str( num_of_drugs                     )  +  "\n"     )
    report_file.write( "    Errors occured: "       + str( error_counter                    )  +  "\n"     )
    report_file.write( "    Percentage of errors: " + str( reliability                      )  +  "\n"     )
    report_file.write( "    Still need to do: "     + str( DRUGS_THAT_OUTPUT_LISTS )  +  "\n\n\n" )


    report_file.close()
    


def ScrapWebpage(file):
    #https://www.ncbi.nlm.nih.gov/pccompound/?term=lovastatin
    #https://www.ncbi.nlm.nih.gov/pcsubstance?term=lovastatin

    # Load in the workbook
    book = xlrd.open_workbook(file)

    #print book.sheet_names() 

    sheet = book.sheet_by_index(CURRENT_SHEET)
    #print sheet.name

    drugs_array = sheet.col_values(COLUMN_CONTAINING_NAMES)
    num_of_drugs = len(drugs_array)

    file_name = sheet.name + "_ids_raw_scrap_sid_full.txt"

    file = open(file_name, "w")
    error_counter = 0.0


    i = 0
    for drug in drugs_array:
        if i != 0:
            try: 
                #compound result if the page exists 
                cid_page = getresult("https://www.ncbi.nlm.nih.gov/pcsubstance/?term=%s" % drug)
                soup = BeautifulSoup(cid_page, 'html.parser')
                desc = soup.findAll(attrs={"name":"description"}) 
                contents = desc[0]['content'].encode('utf-8') 
                meta_list = contents.split()
                #print meta_list
                for idx, ele in enumerate(meta_list):
                    if ele == "SID":
                        cid_actual = meta_list[(idx + 1) % len(meta_list)]
                        print "Drug broad cpd id: ", drug
                        print "SID: ", cid_actual
                        print

                        file.write(drug + ";" + str(cid_actual) + "\n")
            except:
                print "****PAGE DOES NOT EXIST****\n"
                file.write(drug + ";ERROR\n")
                error_counter += 1


                    
        i += 1





    file.close()

    reliability = error_counter / num_of_drugs
    print "Errors occured: ", error_counter
    print "Number of drugs: ", num_of_drugs
    print "Percentage of errors: ", reliability


def clean(data):
    #print type(data)
    if type(data) == float:
        return str(int(data)).strip()
    else:
        return str(data).strip()





if __name__ == "__main__":
    # Quick check before starting
    CheckBeforeRunning( NAME_OF_FILE )
    ListFileSheetColumn()


    #ParseExcel( "Drug_Conversion_test.xlsx" )
    #ScrapWebpage( "Drug_Conversion_test.xlsx" )
    #MakeMaster( "Drug_Conversion_test.xlsx" )



    #testing
    #print StripNames("Afatinib + Cetuximab (ERBB2 inhibitor&EGFR inhibitor 2nd gen + EGFR mAb inhibitor)")
    #print Check4List( NameToCid("romidepsin") )
    #print Check4List( NameToCid("poziotinib") )
