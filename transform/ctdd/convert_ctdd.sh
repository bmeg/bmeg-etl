#!/usr/bin/env bash

# Sample ctdd converter script. Please set-up with your local paths and options.

CWD=$(pwd)
BDIR=$( cd $(dirname $0) && pwd )
CONVERTER_SCRIPT=${BDIR}/convert-ctdd.py

DATADIR=$1
OUTPATH="${CWD}/ctdd.json"

RESPONSE="${DATADIR}/v20.data.curves_post_qc.txt"
METADRUG="${DATADIR}/v20.meta.per_compound.txt"
METACCL="${DATADIR}/v20.meta.per_cell_line.txt"
METAEXPERIMENT="${DATADIR}/v20.meta.per_experiment.txt"
AVEDATA="${DATADIR}/v20.data.per_cpd_avg.txt"

PUBCHEM="${CWD}/ctdd_pubchem.table"

echo $CWD $CONVERTER_SCRIPT $DATADIR $OUTPATH

python $CONVERTER_SCRIPT \
--response $RESPONSE \
--metadrug $METADRUG \
--metacellline $METACCL \
--metaexperiment $METAEXPERIMENT \
--pubchem $PUBCHEM \
--data $AVEDATA \
--multi $OUTPATH \
--format json
