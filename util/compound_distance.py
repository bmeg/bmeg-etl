#!/usr/bin/env python

import gzip
import numpy as np
import json
from multiprocessing import Pool
import sys
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

NCPUS=8

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
    return None, None

def pairList(compounds, fingerprints):
    for i in compounds:
        if i in fingerprints:
            c1 = fingerprints[i]
            for j in compounds:
                if j in fingerprints:
                    c2 = fingerprints[j]
                    yield i, c1, j, c2

def calcDistance(args):
    i, c1, j, c2 = args
    tan_sim = DataStructs.FingerprintSimilarity(c1, c2, metric=DataStructs.TanimotoSimilarity)
    dice_sim = DataStructs.FingerprintSimilarity(c1, c2, metric=DataStructs.DiceSimilarity)
    if dice_sim >= 0.5:
        return json.dumps( {
            "compound_1" : i,
            "compound_2" : j,
            "morgan_fingerprint_2_dice":dice_sim,
            "morgan_fingerprint_2_tanimoto":tan_sim
        } )
    return None


if __name__ == "__main__":

    records={}
    fingerprints = {}
    with Pool(processes=NCPUS) as pool:
        for id, fingerprint in pool.imap_unordered(lineProcess, sys.stdin):
            if id is not None:
                fingerprints[id] = fingerprint

    # Calculate similarity
    compounds = list(fingerprints.keys())

    with Pool(processes=NCPUS) as pool:
        for line in pool.imap_unordered(calcDistance, pairList(compounds, fingerprints)):
            if line is not None:
                print(line)