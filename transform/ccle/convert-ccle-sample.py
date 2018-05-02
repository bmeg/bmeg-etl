#!/usr/bin/env python

import json
import sys
import csv
from bmeg import clinical_pb2
# from ga4gh.bio_metadata_pb2 import Biosample
from google.protobuf import json_format

"""

Clinical data download
curl -o CCLE_sample_info_file_2012-10-18.txt "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_22/CCLE_sample_info_file_2012-10-18.txt?downloadff=true&fileId=6801"

Columns:
    CCLE name
    Cell line primary name
    Cell line aliases
    Gender
    Site Primary
    Histology
    Hist Subtype1
    Notes
    Source
    Expression arrays
    SNP arrays
    Oncomap
    Hybrid Capture Sequencing

Example:
    1321N1_CENTRAL_NERVOUS_SYSTEM
    1321N1

    M
    central_nervous_system
    glioma
    astrocytoma
    Identical lines: U-118 MG, U-138 MG and 1321N1 share high SNP identity
    ECACC
    NIECE_p_NCLE_RNA3_HG-U133_Plus_2_B06_296024
    HONEY_p_NCLE_DNAAffy3_S_GenomeWideSNP_6_E09_293392
    yes

"""

with open(sys.argv[1]) as handle, open(sys.argv[2], "w") as out_handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for line in reader:
        msg = clinical_pb2.Biosample()
        msg.id = line['CCLE name']
        msg.dataset_id = 'ccle'
        msg.name = line['CCLE name']
        msg.source = 'ccle'

        # Does not exist? msg.individual_id =
        # Does not exist? msg.created = 
        # Does not exist? msg.updated =

        for k, v in line.items():
            msg.attributes[k] = v
        #print json.dumps(line)
        d = json_format.MessageToDict(msg)
        out_handle.write(json.dumps(d) + "\n")
