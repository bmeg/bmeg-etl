#!/usr/bin/env python

import polars as pl
import pandas as pd
import argparse
import glob
import string
import json
import sys
## This program is written specifically for matching the drug names in the data emitted from the transform/pharmacodb/drug_response_tsv.R file and the tranform/g2p/transform.yaml from June 2023 to CHEMBL IDs. Currently it requires three input files: the drug response tsv, the treatment table with some synonyms, and the CHEMBL id-synonym table produced by the chembl synonyms.yaml file. Dependencies include polars and pandas python libraries. WHEN MODIFYING THIS PROGRAM, BE AWARE THAT TWO SEPARATE PARTS OF BMEG USE IT

parser = argparse.ArgumentParser(description="Add chembl id's to a tsv or json containing treatment names")
parser.add_argument("-r", "--response", help="Drug response/phenotype file. Use ... to represent wildcards")
parser.add_argument("-c", "--chembl", help="Chembl matchup table")
parser.add_argument("-i", "--intermediate", help="Intermediate table file")
args = parser.parse_args()

gpath = args.response.replace('...', '*')

responseFiles = glob.glob(gpath)

def matchChembl(treatmentFile, chemblFile, intermediateFile):
  #parser.add_argument("-rd", "--response_directory", help="Drug response/phenotype")
  mapping = {ord(c): None for c in string.whitespace}
  
  #importing main tsv:
  if treatmentFile[-7:] == ".tsv.gz" or treatmentFile[-4:] == '.tsv':
    pharmaco = pl.read_csv(treatmentFile, separator="\t", infer_schema_length=0)
  elif treatmentFile[-7:] == "json.gz":
    pharmaco = pl.from_pandas(pd.read_json(treatmentFile, lines=True))
  else:
    exit("Unsupported file extension " + treatmentFile[-7:] + " drug response can only be .tsv or .json")
  pharmacofinal = pharmaco.select(pl.exclude("UNIQUEtreatmentid")).with_columns(pl.col("treatmentid").apply(lambda x: x.lower().translate(mapping)).alias("standard"))
  
  #importing pharmaco curation table
  pharmacoTreat = pl.from_pandas(pd.read_json(intermediateFile, lines=True)) 
  treatfinal = pharmacoTreat.with_columns(pl.col("PROJECTtreatmentid").apply(lambda x: x.lower().translate(mapping)).alias("standard")).unique("standard") 
  
  #importing reference for chembl ids
  chembl = pl.from_pandas(pd.read_json(chemblFile, lines=True))
  chemblfinal = chembl.with_columns(pl.col("synonyms").apply(lambda x: x.lower().translate(mapping)).alias("standard")).rename({"id":"chID"}).sort("chID").unique("standard", keep="first")
  
  combine = pharmacofinal.join(treatfinal, on="standard", how="left")
  
  combinefinal = combine.with_columns(pl.coalesce("UNIQUEtreatmentid","treatmentid").apply(lambda x: x.lower().translate(mapping)).alias("standard"))
  
  ct = combinefinal.join(chemblfinal, on="standard", how="left")
  #.unique(["responseID"])
  
  ctfinal = ct.with_columns(pl.coalesce(["chID","standard"]).apply(lambda x: x.upper() if x.upper()[0:6] == "CHEMBL" else None).alias("chemblid"))
  
  for row in ctfinal.select(pl.exclude('compounds')).iter_rows(named=True):
    rowd = {k:v for k,v in row.items() if v is not None}
    #if "sampleid" in rowd and "treatmentid" in rowd:
    print(json.dumps(rowd))

for rFile in responseFiles:
  matchChembl(rFile, args.chembl, args.intermediate)
