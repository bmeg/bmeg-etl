#! /usr/bin/python
'''
Authors: Malisa Smith smimal@ohsu.edu, Ryan Spangler spanglry@ohsu.edu

This program converts CTDD drug response information into
protobuf data based on the BMEG sample.proto schema.

Source: ftp://caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip

The four files of interest (for this converter) are:
1) v20.data.curves_post_qc.txt
2) v20.meta.per_compound.txt
3) v20.meta.per_cell_line.txt
4) v20.meta.per_experiment.txt
'''

import tempfile
import subprocess
from google.protobuf import json_format
import json, sys, argparse, os
import csv #for drug data
import string
import re
import pandas
from bmeg import phenotype_pb2, conversions

def parse_args(args):
    # We don't need the first argument, which is the program name
    args = args[1:]

    # Construct the parser
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Now add all the options to it
    parser.add_argument('--zip', type=str, help='Path to zip archive', default=None)
    parser.add_argument('--response', type=str, help='Path to the drug response experiment data you want to import')
    parser.add_argument('--metadrug', type=str, help='Path to the drug meta data you want to import')
    parser.add_argument('--metacellline', type=str, help='Path to the cell line meta data you want to import')
    parser.add_argument('--metaexperiment', type=str, help='Path to the experiment meta data you want to import')
    parser.add_argument('--pubchem', type=str, help='Path to file to map drug names to pubchem ids')
    parser.add_argument('--data', type=str, help='Path to the experiment data you want to import')
    parser.add_argument('--out', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--multi', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')
    args = parser.parse_args(args)
    
    if args.zip is not None:
        workdir = tempfile.mkdtemp(prefix="ctrp-data", dir="./")
        subprocess.check_call("unzip %s" % (os.path.abspath(args.zip)), cwd=workdir, shell=True)
        args.response = os.path.join(workdir, "v20.data.curves_post_qc.txt")
        args.metadrug = os.path.join(workdir, "v20.meta.per_compound.txt")
        args.metacellline = os.path.join(workdir, "v20.meta.per_cell_line.txt")
        args.metaexperiment = os.path.join(workdir, "v20.meta.per_experiment.txt")
        args.data = os.path.join(workdir, "v20.data.per_cpd_avg.txt")
    
    return args

########################################

def find_biosample(state, source, barcode, sample_type):
    sample_name = 'ccle:' + barcode
    biosample = state['Biosample'].get(sample_name)
    if biosample is None:
        biosample = schema.Biosample()
        biosample.name = sample_name
        biosample.dataset_id = "CCLE"
        biosample.source = source
        biosample.barcode = barcode
        biosample.sampleType = sample_type
        state['Biosample'][sample_name] = biosample

    return biosample

def append_unique(l, i):
    if not i in l:
        l.append(i)

def get_compounds(name, name_table):
    res = re.search(r'(.*):(.*) \((\d+):(\d+)', name)
    if res:
        groups = res.groups()
        name1 = conversions.pubchem(groups[0])
        name2 = conversions.pubchem(groups[1])

        # if name_table[groups[0]] != "none":
        #     name1 = name_table[groups[0]]
        # else:
        #     name1 = groups[0]
        # if name_table[groups[1]] != "none":
        #     name2 = name_table[groups[1]]
        # else:
        #     name2 = groups[1]

        return [(name1, groups[2]), (name2, groups[3])]
    else:
        return [(conversions.pubchem(name), 1)]

        # if name in name_table and name_table[name] != 'none':
        #     return [(name_table[name], 1)]
        # else:
        #     return [(name, 1)]

def process_drugs(emit, input, name_table): #row is a namedtuple
    compounds = set()
    
    for row in input.itertuples():
        # create drug message for CTDD compound
        for compound_name, ratio in get_compounds(row.cpd_name, name_table):
            if compound_name not in compounds:
                compound = phenotype_pb2.Compound()
                compound.id = compound_name
                # compound.gid = compound_name
                compound.name = row.cpd_name
                compound.smiles = row.cpd_smiles
                compound.status = row.cpd_status
                target = row.gene_symbol_of_protein_target
                if target and isinstance(target, str):
                    targets = target.split(';')
                    for hugo in targets:
                        ensembl = conversions.hugo_ensembl(hugo)
                        if ensembl != '':
                            compound.target.append(ensembl)

                compound.report = row.target_or_activity_of_compound
                compound.rationale = row.inclusion_rationale
                compound.synonyms.append(row.broad_cpd_id)
                emit(compound)
                compounds.add(compound_name)

def process_response(emit, input, data, name_table):
    gid_set = set()
    for row in input.itertuples():
        if isinstance(row.ccle_primary_site, str):
            site = row.ccle_primary_site.upper()
            sample_name = '%s_%s' % (row.ccl_name, site)
            sample = sample_name
            gid = "ccle-response:%s/%s" % (sample_name, row.cpd_name)

            if gid not in gid_set:
                response = phenotype_pb2.ResponseCurve()
                # response.gid = gid
                response.source = 'ccle'
                response.responseType = phenotype_pb2.ResponseCurve.ACTIVITY
                compound_set = get_compounds(row.cpd_name, name_table)
                compound_total = sum( float(a[1]) for a in compound_set )
                for compound_name, ratio in compound_set:
                    n = response.compounds.add()
                    n.compound = compound_name
                    n.ratio = float(ratio) / compound_total
                response.sample = sample

                s = response.summary.add()
                s.type = phenotype_pb2.ResponseSummary.EC50
                s.value = row.apparent_ec50_umol
                s.unit = "uM"
            
                s = response.summary.add()
                s.type = phenotype_pb2.ResponseSummary.AUC
                s.value = row.area_under_curve
                s.unit = "uM"

                for m in data.loc[lambda x: x.master_cpd_id==row.master_cpd_id, : ].loc[lambda x: x.experiment_id==row.experiment_id].itertuples():
                    dr = response.values.add()
                    dr.dose = m.cpd_conc_umol
                    dr.response = m.cpd_expt_avg_log2
            
                emit(response)
                gid_set.add(gid)

def convert_all_ctdd(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, pubchemPath, dataPath, out, multi=None):
    # Read in Compound information into a pandas dataframe.
    compound_df = pandas.read_table(metadrugPath)
    # Read in Cell line information
    ccl_df = pandas.read_table(metacelllinePath)
    # Read in data curves for experiments
    datacurves_df = pandas.read_table(responsePath)
    # Read in meta experimental data
    metaexperiment_df = pandas.read_table(metaexperimentPath)
    
    ctdd_merged = pandas.merge(datacurves_df, metaexperiment_df, how='left', on=['experiment_id']) # merge experiment data
    ctdd_merged = pandas.merge(ctdd_merged, compound_df, how='left', on=['master_cpd_id']) # merge with compound data frame
    ctdd_merged = pandas.merge(ctdd_merged, ccl_df, how='left', on=['master_ccl_id']) # merge with cell line data frame
    
    ctdd_data = pandas.read_table(dataPath)
    #print ctdd_merged
    
    out_handles = {}
    def emit_json_single(message):
        if 'main' not in out_handles:
            out_handles['main'] = open(out, "w")
        msg = json.loads(json_format.MessageToJson(message))
        msg["#label"] = message.DESCRIPTOR.full_name
        out_handles['main'].write(json.dumps(msg))
        out_handles['main'].write("\n")

    def emit_json_multi(message):
        if message.DESCRIPTOR.full_name not in out_handles:
            out_handles[message.DESCRIPTOR.full_name] = open(multi + "." + message.DESCRIPTOR.full_name + ".json", "w")
        msg = json.loads(json_format.MessageToJson(message))
        out_handles[message.DESCRIPTOR.full_name].write(json.dumps(msg))
        out_handles[message.DESCRIPTOR.full_name].write("\n")

    if out is not None:
        emit = emit_json_single
    if multi is not None:
        emit = emit_json_multi
    
    ctdd_merged.to_csv("test.out", sep="\t")    
    
    pubchem_ids = {}
    with open(pubchemPath) as handle:
        for line in handle:
            row = line.rstrip().split("\t")
            if row[1] != 'none':
                pubchem_ids[row[0]] = row[1] # "pubchem:" + row[1]
            else:
                pubchem_ids[row[0]] = row[1] # "compound" + row[1]
    process_drugs(emit, ctdd_merged, pubchem_ids)
    process_response(emit, ctdd_merged, ctdd_data, pubchem_ids)
    
########################################

def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    msg['#label'] = message.DESCRIPTOR.name
    return json.dumps(msg)


def convert_to_profobuf(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, pubchemPath, dataPath, out, multi):
    if responsePath and metadrugPath and metacelllinePath and metaexperimentPath and (out or multi) and format:
        convert_all_ctdd(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, pubchemPath, dataPath, out, multi)
    else:
        print("Please include all arguments")

    #write_messages(state, outpath, format)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_to_profobuf(options.response, options.metadrug, options.metacellline, options.metaexperiment, options.pubchem, dataPath=options.data, out=options.out, multi=options.multi)
