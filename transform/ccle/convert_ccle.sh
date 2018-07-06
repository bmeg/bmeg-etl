#!/usr/bin/env bash

# Sample ccle converter script. Please set-up with your local paths and options.

CWD=$(pwd)
CONVERTER_SCRIPT=${CWD}/convert-ccle.py

OUTPATH="${CWD}/tmp_converted/ccle_converted.json"

MAF="${CWD}/raw-data/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf"
CSV="${CWD}/raw-data/CCLE_NP24.2009_Drug_data_2015.02.24.csv"

python $CONVERTER_SCRIPT --maf $MAF --csv $CSV --out $OUTPATH --format json
