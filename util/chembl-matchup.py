#!/usr/bin/env python

import polars as pl
import pandas as pd
import argparse
import string
import json

# This program is written specifically for matching the drug names in the data emitted from the 
# transform/pharmacodb/drug_response_tsv.R file and the tranform/g2p/transform.yaml from 
# June 2023 to CHEMBL IDs. Currently it requires three input files: the drug response tsv, the 
# treatment table with some synonyms, and the CHEMBL id-synonym table produced by the chembl synonyms.yaml file. 
# Dependencies include polars and pandas python libraries. 
# WHEN MODIFYING THIS PROGRAM, BE AWARE THAT TWO SEPARATE PARTS OF BMEG USE IT


def matchChembl(responseFile, chemblFile, treatmentFile):
  mapping = {ord(c): None for c in string.whitespace}
  
  #importing main tsv:
  if responseFile.endswith(".tsv.gz") or responseFile.endswith('.tsv'):
    pharmaco = pl.read_csv(responseFile, separator="\t", infer_schema_length=0)
  elif responseFile.endwith("json.gz"):
    pharmaco = pl.from_pandas(pd.read_json(responseFile, lines=True))
  else:
    exit("Unsupported file extension " + responseFile + " drug response can only be .tsv or .json")
  pharmacofinal = pharmaco.select(pl.exclude("UNIQUEtreatmentid")).with_columns(pl.col("treatmentid").apply(lambda x: x.lower().translate(mapping)).alias("standard"))
  
  #importing pharmaco curation table
  pharmacoTreat = pl.from_pandas(pd.read_json(treatmentFile, lines=True)) 
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

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Add chembl id's to a tsv or json containing treatment names")
  parser.add_argument("-c", "--chembl", help="Chembl matchup table")
  parser.add_argument("-t", "--treatment", help="Treatment table file")
  parser.add_argument("response", nargs="+", help="Drug response/phenotype file")
  args = parser.parse_args()

  for res in args.response:
    matchChembl(res, args.chembl, args.treatment)
