#!/usr/bin/env python

import gzip
import numpy as np
import json
import functools
from multiprocessing import Pool
import sys
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import argparse


fingerprints = {}
compounds = []

def init_worker(fingerprintsIn, compoundsIn):
    global fingerprints
    global compounds
    fingerprints = fingerprintsIn
    compounds = compoundsIn

def list2fingerprint(fingerprint):
    '''input fingerprint (stored as list in record) and output morgan fingerprint from rdskit
    '''
    fp = np.array(fingerprint)
    bitstring="".join(fp.astype(str))
    fp_bit = DataStructs.cDataStructs.CreateFromBitString(bitstring)
    return fp_bit

def lineProcess(line):
    d = json.loads(line)
    if "morgan_fingerprint_2" in d:
        return d["chembl_id"], list2fingerprint(d["morgan_fingerprint_2"])
        #fp = np.array(d["morgan_fingerprint_2"])
        #bitstring="".join(fp.astype(str))
        #return d["chembl_id"], bitstring
    return None, None

def lineProcessFHIR(line):
    d = json.loads(line)
    if "classification" in d.keys() and d["classification"]:
        morgan_fingerprint = json.loads(d['classification'][0]['coding'][0]['display'])
        return d["id"], list2fingerprint(morgan_fingerprint)
    return None, None

def calcDistance(n):
    out = []
    i = compounds[n]
    c1 = fingerprints[i]
    #c1 = DataStructs.cDataStructs.CreateFromBitString(fingerprints[i])
    for j in compounds[n+1:]:
        if j in fingerprints:
            c2 = fingerprints[j]
            #c2 = DataStructs.cDataStructs.CreateFromBitString(fingerprints[j])
            tan_sim = DataStructs.FingerprintSimilarity(c1, c2, metric=DataStructs.TanimotoSimilarity)
            dice_sim = DataStructs.FingerprintSimilarity(c1, c2, metric=DataStructs.DiceSimilarity)
            if dice_sim >= 0.7:
                o = json.dumps( {
                    "compound_1" : i,
                    "compound_2" : j,
                    "morgan_fingerprint_2_dice":dice_sim,
                    "morgan_fingerprint_2_tanimoto":tan_sim
                } )
                out.append(o)
    return out


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--ncpus", default=8, type=int)
    parser.add_argument("-s", "--shard", default=1, type=int)
    parser.add_argument("-t", "--total-shards", default=1, type=int)
    parser.add_argument("-i", "--input", default=None)
    parser.add_argument("-o", "--output", default=None)

    args = parser.parse_args()

    input = sys.stdin
    if args.input is not None:
        if args.input.endswith(".gz"):
            input = gzip.open(args.input)
        else:
            input = open(args.input)

    with Pool(processes=args.ncpus) as pool:
        for id, fingerprint in pool.imap(lineProcess, input):
            if id is not None:
                fingerprints[id] = fingerprint

    for i in fingerprints.keys():
        compounds.append(i)

    output = sys.stdout
    if args.output is not None:
        if args.output.endswith(".gz"):
            output = gzip.open(args.output, "wt")
        else:
            output = open(args.output, "wt")
    with Pool(initializer=init_worker,  initargs=[fingerprints, compounds], processes=args.ncpus) as pool:
        for lines in pool.imap(calcDistance, list( i+args.shard-1 for i in  range(0, len(compounds)-(args.shard-1), args.total_shards))):
            if len(lines):
                for l in lines:
                    output.write(l)
                    output.write("\n")
    output.close()
