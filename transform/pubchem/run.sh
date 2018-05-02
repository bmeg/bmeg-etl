#!/usr/bin/env bash


PREFIX="ctdd.json.bmeg.Compound.json"
RCPREFIX="ctdd.json.bmeg.ResponseCurve.json"

INPUTPATH="protograph_files/pre_Protograph"
OUTPUTPATH="protograph_files/post_Protograph"

echo "Starting Protobuf prepping..."
#./convert_ctdd.sh
echo "Moving files into approapriate dir..."
#mv ctdd.json.bmeg.Compound.json $INPUTPATH
echo "  File Moved!"
echo ""
echo "Processing protobuf output into protograph input..."
#java -jar protograph.jar --protograph protograph.yml \
#--input $INPUTPATH/$PREFIX \
#--output $OUTPUTPATH/$PREFIX
java -jar protograph.jar --protograph protograph.yml \
--input $INPUTPATH/$RCPREFIX \
--output $OUTPUTPATH/$RCPREFIX

echo ""
echo "Finished processing! Ready to go for ingestion!"
echo ""
echo "Have a great day!"
