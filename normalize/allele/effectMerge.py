#!/usr/bin/env python

from bx import (intervals, misc)
import pandas as pd
import polars as pl
import argparse
import json

k = 100
parser = argparse.ArgumentParser(description="Basically annotates a vcf of allele effects using the amino acid equivalent of a bed file")
parser.add_argument("-f", "--featuresFile", help="json.gz of features to match up (has transcript, start and end cols)")
parser.add_argument("-v", "--variantsFile", help="json.gz of variants to match up (has ensembl_transcript and aa_position cols")
args = parser.parse_args()

featuresTotal = pl.from_pandas(pd.read_json(args.featuresFile, lines=True)).unique()
variants = pl.from_pandas(pd.read_json(args.variantsFile, lines=True)).unique()
var1Total = variants.filter(pl.col('aa_ref').is_not_null())
transcripts = var1Total['ensembl_transcript'].unique()
features = featuresTotal.filter(pl.col('transcript').is_in(transcripts).and_(pl.col('start').is_not_null(),pl.col('end').is_not_null()))
var1 = var1Total

intersecters = {}

for row in features.iter_rows(named=True):
  if row['transcript'] not in intersecters:
    intersecters[row['transcript']] = intervals.Intersecter()
  intersecters[row['transcript']].add_interval(intervals.Interval(int(row['start']), int(row['end'])))

rowFeatures = []

for row in var1.iter_rows(named=True):
  if row['ensembl_transcript'] in intersecters:
    rowFeatures.append({'aa_position': row['aa_position'], 'transcript': row['ensembl_transcript'], 'intervals': str(intersecters[row['ensembl_transcript']].find(row['aa_position'], row['aa_position']))[1:-1]})

intermediate = pl.from_dicts(rowFeatures).unique().filter(pl.col('intervals').ne(''))
iFinal = intermediate.with_columns(pl.col('intervals').str.slice(8).str.split(', Interval').alias('range')).explode('range').drop('intervals')
featuresFinal = features.with_columns(('('+pl.col('start').cast(int).cast(str)+ ', '+pl.col('end').cast(int).cast(str)+')').alias('range')).drop(['start', 'end'])

iTest = iFinal.join(featuresFinal, on=['transcript','range'], how='left').drop('range').unique()

preFinal = iTest.groupby(['aa_position', 'transcript']).agg(pl.all().drop_nulls().first())

final = variants.join(preFinal, left_on=['aa_position', 'ensembl_transcript'],right_on=['aa_position','transcript'], how='left')




#while len(var1) > 0:
#  merge1 = var1[:k].join(features, left_on='ensembl_transcript', right_on='transcript', how='left')
#  var1 = var1[k:]
#  merge2 = merge1.filter(pl.col('aa_position').is_between(pl.col('start'),pl.col('end'))).sort('id').select(pl.all().exclude('start','end'))
#  
#  merge3 = merge2.groupby('id').agg(pl.all().drop_nulls().unique().first())
#  #merge2.groupby('id').agg(pl.all().drop_nulls().explode().unique().cast(pl.Utf8).str.concat(', '))
#  
#  for row in merge3.iter_rows(named=True):
#    rowd = {k:v for k,v in row.items() if v is not None}
#    #if "sampleid" in rowd and "treatmentid" in rowd:
#    print(json.dumps(rowd))
#
for row in final.iter_rows(named=True):
  rowd = {k:v for k,v in row.items() if v is not None}
  #if "sampleid" in rowd and "treatmentid" in rowd:
  print(json.dumps(rowd))
