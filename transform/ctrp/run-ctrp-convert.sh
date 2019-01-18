#!/bin/bash

TMP=`mktemp -d -p .`
unzip -d $TMP source/ctrp/CTRPv2.0_2015_ctd2_ExpandedDataset.zip
python transform/ctrp/ctrp-convert.py $TMP
rm -rf $TMP
