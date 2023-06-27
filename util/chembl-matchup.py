#!/usr/bin/env python

import polars as pl
import pandas as pd
import string
import json
import sys
## This program is written specifically for matching the drug names in the data emitted from the drug_response_tsv.R file from June 2023 to CHEMBL IDs. Currently it requires three input files: the drug response tsv, the treatment table with some synonyms, and the CHEMBL id-synonym table produced by the chembl synonyms.yaml file.

mapping = {ord(c): None for c in string.whitespace}

#importing main tsv:


if sys.argv[1][-4:] == ".tsv":
  pharmaco = pl.read_csv(sys.argv[1], separator="\t", infer_schema_length=0)
elif sys.argv[1][-4:] == "json":
  pharmaco = pl.from_pandas(pd.read_json(sys.argv[1], lines=True))
else:
  exit("Unsupported file extension " + sys.argv[1][-4:] + " drug response can only be .tsv or .json")
pharmacofinal = pharmaco.select(pl.exclude("UNIQUEtreatmentid")).with_columns(pl.col("treatmentid").apply(lambda x: x.lower().translate(mapping)).alias("standard"))

#importing pharmaco curation table
pharmacoTreat = pl.read_csv(sys.argv[2], separator="\t", infer_schema_length=0)
treatfinal = pharmacoTreat.with_columns(pl.col("PROJECTtreatmentid").apply(lambda x: x.lower().translate(mapping)).alias("standard")).unique() 

#importing reference for chembl ids
chembl = pl.read_ndjson(sys.argv[3])
chemblfinal = chembl.with_columns(pl.col("synonyms").apply(lambda x: x.lower().translate(mapping)).alias("standard"))

combine = pharmacofinal.join(treatfinal, on="standard", how="left").unique("responseID")


combinefinal = combine.with_columns(pl.coalesce("UNIQUEtreatmentid","treatmentid").apply(lambda x: x.lower().translate(mapping)).alias("standard"))

ct = combinefinal.join(chemblfinal, on="standard", how="left").unique(["responseID"])

ctfinal = ct.with_columns(pl.coalesce(["id","standard"]).apply(lambda x: x.upper() if x.upper()[0:6] == "CHEMBL" else None).alias("chemblid"))

for row in ctfinal.iter_rows(named=True):
  rowd = {k:v for k,v in row.items() if v is not None} 
  if "sampleid" in rowd and "treatmentid" in rowd:
    print(json.dumps(rowd))

