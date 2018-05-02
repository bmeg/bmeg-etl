#!/usr/bin/env python

import re
import argparse
import json
import gzip
from phenotype_pb2 import Compound, Assay, ResponseCurve
from urllib2 import urlopen
from urllib import quote
from google.protobuf import json_format
import csv
import sys


#from zeep import Client

def message_to_json(message, indent=None):
    if indent is None:
        return json.dumps(json_format.MessageToDict(message))
    return json.dumps(json_format.MessageToDict(message), indent=indent)

RECORD = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/record/JSON"
ASSAY_SUMMARY = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON"
ASSAYS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/aids/JSON"
SYNONYMS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/JSON"
DESCRIPTION = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/description/JSON"
PUBMED = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/xrefs/pubmedid/JSON"
NAME_CID = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/JSON"

SIDS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/sids/JSON"

ASSAY_RECORD = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/record/JSON"
DOSE_RESPONSE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/doseresponse/JSON"

SID2CID = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/%s/cids/JSON"


def run_search(args):
    u = NAME_CID % (quote(args.name))
    try:
        record_txt = urlopen(u).read()
    except:
        return
    record = json.loads(record_txt)
    record['Name'] = args.name
    print json.dumps(record)

def run_compound_extract(args):
    """
    #ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
    GENE_FILE = sys.argv[1]
    gene_map = {}
    with gzip.GzipFile(GENE_FILE) as handle:
        for line in handle:
            row = line.rstrip().split("\t")
            gene_map[row[1]] = row[2]

    CID = sys.argv[2]

    record_txt = urlopen(RECORD % CID).read()
    print record_txt
    """

    compound_handle = open(args.out, "a+")

    for CID in args.ids:
        out = Compound()
        out.id = 'CID{}'.format(CID)

        record_txt = urlopen(RECORD % CID).read()
        record = json.loads(record_txt)
        for compound in record["PC_Compounds"]:
            #print compound.keys()
            props = {}
            for prop in compound["props"]:
                if "name" in prop["urn"]:
                    label = prop["urn"]["label"] + ":" +  prop["urn"]["name"]
                else:
                    label = prop["urn"]["label"]
                for i in prop["value"].values():
                    props[label] = i
            #print props
        #print record
        out.smiles = props["SMILES:Canonical"]

        syn_txt = urlopen(SYNONYMS % CID).read()
        syn = json.loads(syn_txt)

        for info in syn["InformationList"]["Information"]:
            for name in info["Synonym"]:
                if name.startswith("CHEBI:"):
                    out.chebi_id = name
                out.synonyms.append(name)

        desc_txt = urlopen(DESCRIPTION % CID).read()
        desc = json.loads(desc_txt)
        for info in desc["InformationList"]["Information"]:
            if 'Title' in info:
                out.name = info['Title']
            if 'Description' in info:
                out.description.append(info['Description'])
        #print desc

        xref_txt = urlopen(PUBMED % CID).read()
        xref = json.loads(xref_txt)
        #print xref
        for info in xref["InformationList"]["Information"]:
            if 'PubMedID' in info:
                for pmid in info['PubMedID']:
                    out.pubmed.append(str(pmid))

        aids_txt = urlopen(ASSAYS % CID).read()
        aids = json.loads(aids_txt)
        for info in aids["InformationList"]["Information"]:
            if 'AID' in info:
                for aid in info['AID']:
                    out.assays.append(str(aid))

        sids_txt = urlopen(SIDS % CID).read()
        sids = json.loads(sids_txt)
        for info in aids["InformationList"]["Information"]:
            if 'SID' in info:
                for sid in info['SID']:
                    out.sids.append(str(sid))

        compound_handle.write( message_to_json(out) )
        compound_handle.write( "\n" )


        """
        assay_txt = None
        try:
            assay_txt = urlopen(ASSAY_SUMMARY % CID).read()
        except:
            pass

        if assay_txt is not None:
            assay = json.loads(assay_txt)

            columns = assay["Table"]["Columns"]["Column"]

            for row in assay["Table"]["Row"]:
                cells = row["Cell"]
                data = dict(zip(columns, cells))
                if len(data['Target GeneID']):
                    print data
        """
    compound_handle.close()
    #chebi_id = sys.argv[1]
    #client = Client('http://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl')
    #print client.service.getCompleteEntity(chebi_id)

def run_assay_extract(args):

    for aid in args.ids:
        aidrecord_txt = None
        try:
            aidrecord_txt = urlopen(ASSAY_RECORD % aid).read()
        except:
            pass
        if aidrecord_txt is not None:
            aidrecord = json.loads(aidrecord_txt)
            #print aidrecord_txt

            out = Assay()
            for record in aidrecord['PC_AssayContainer']:
                out.id = str(record['assay']['descr']['aid']['id'])
                out.name = record['assay']['descr']['name']

                for cmt in record['assay']['descr']['comment']:
                    res = re.search(r'Cell Line:\s*(\w+)', cmt)
                    if res:
                        out.cellline = res.group(1)

                for xref in record['assay']['descr']["xref"]:
                    if 'pmid' in xref['xref']:
                        out.pubmed.append( str(xref['xref']['pmid']) )

            dose_txt = None
            try:
                dose_txt = urlopen(DOSE_RESPONSE % aid).read()
            except:
                pass
            if dose_txt is not None:
                dose = json.loads(dose_txt)
                columns = dose["Table"]["Columns"]["Column"]
                dose_out = {}
                for line in dose["Table"]["Row"]:
                    row = dict(zip(columns, line["Cell"]))
                    sid = row["SID"]
                    if sid not in dose_out:
                        dose_out[sid] = out.drugresponse.add()

                        try:
                            cid_txt = urlopen(SID2CID % (sid)).read()
                            cid_res = json.loads(cid_txt)
                            #print "howdy", cid_res['InformationList']['Information']
                            #print cid_res['InformationList']['Information'][0]['CID'][0]
                            cid = str(cid_res['InformationList']['Information'][0]['CID'][0])
                            c = dose_out[sid].compounds.add()
                            c.compound = "CID:%s" % (cid)
                        except:
                            c = dose_out[sid].compounds.add()
                            c.compound = "SID:%s" % (sid)

                    o = dose_out[sid].values.add()
                    o.dose = float(row["Concentration"])
                    o.response = float(row["Response"])
                    #print row


            print message_to_json(out)


def run_table_xform(args):
    """ Extract ids from mapping tables.
        Redirect to run_compound_extract()
        Misses create """
    TABLES = [
        'https://raw.githubusercontent.com/biostream/gdc-transform/master/tcga_pubchem.map',
        'https://raw.githubusercontent.com/biostream/ctdd-transform/master/ctdd_pubchem.table',
        'https://raw.githubusercontent.com/biostream/ccle-transform/master/ccle_pubchem.txt',
        'https://raw.githubusercontent.com/biostream/gdsc-transform/master/gdsc_pubchem.table',
    ]
    for url in TABLES:
        sys.stderr.write('run_table_xform processing {}\n'.format(url))
        tsvin = urlopen(url)
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            name = row[0]
            id = row[1]
            # deep copy args
            args_copy = argparse.Namespace(**vars(args))
            args_copy.ids = [id]
            args_copy.out = args.table_xform_out
            try:
                run_compound_extract(args_copy)
            except Exception as e:
                sys.stderr.write('  error {} {} {}\n'
                                 .format(name, id, e))
                # create unknown
                out = Compound()
                out.id = 'UNKNOWN:{}'.format(name)
                out.name = name
                compound_handle = open(args_copy.out, "a+")
                compound_handle.write(message_to_json(out))
                compound_handle.write("\n")
                compound_handle.close()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers()

    parser_search = subparser.add_parser("search")
    parser_search.add_argument("name")
    parser_search.set_defaults(func=run_search)

    parser_compound = subparser.add_parser("compound-extract")
    parser_compound.add_argument("--out", default="compound.json")
    parser_compound.add_argument("ids", nargs="+")
    parser_compound.set_defaults(func=run_compound_extract)

    parser_assay = subparser.add_parser("assay-extract")
    parser_assay.add_argument("--base", default="out")
    parser_assay.add_argument("ids", nargs="+")
    parser_assay.set_defaults(func=run_assay_extract)

    table_xform = subparser.add_parser("table-xform")
    table_xform.add_argument("--table_xform_out",
                             help="json for table xform output",
                             default="compound.json")
    table_xform.set_defaults(func=run_table_xform)

    args = parser.parse_args()
    args.func(args)
