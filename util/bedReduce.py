#!/usr/bin/env python

import polars as pl
import sys
#Simplifying bed files so that the snpEff tool can grab the most relevant info. snpEff only uses the first 4 columns
inpath = sys.argv[1]
outName = '.'.join(inpath.split('.')[0:-1]) + '.4cols.bed'

domains = pl.read_csv(inpath, separator='\t', infer_schema_length=0)
domains = domains.with_columns(pl.coalesce('comments','longName','name').str.replace_all('"','-double-prime').str.replace_all("'", '-prime').alias('info')).select('#chrom','chromStart','chromEnd','info')
domains.write_csv(outName, separator='\t')
