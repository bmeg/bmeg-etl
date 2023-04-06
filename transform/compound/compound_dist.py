#!/usr/bin/env python

import gzip
import numpy as np
import json
import sys
from tqdm import tqdm
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

def list2fingerprint(fingerprint):
    '''input fingerprint (stored as list in BMEG) and output morgan fingerprint from rdskit
    '''
    fp = np.array(fingerprint)
    bitstring="".join(fp.astype(str))
    fp_bit = DataStructs.cDataStructs.CreateFromBitString(bitstring)
    return fp_bit

if __name__ == "__main__":
    inPath = sys.argv[1]

    records={}
    fingerprints = {}
    with gzip.open(inPath) as handle:
        for line in tqdm(handle):
            d = json.loads(line)
            records[d["chembl_id"]] = d
            if "morgan_fingerprint_2" in d:
                fingerprints[d["chembl_id"]] = list2fingerprint(d["morgan_fingerprint_2"])

    # Calculate similarity
    compounds = list(fingerprints.keys())

    for i in compounds:
        if i in fingerprints:
            c1 = fingerprints[i]
            for j in compounds:
                if j in fingerprints:
                    c2 = fingerprints[j]
                    tan_sim = DataStructs.FingerprintSimilarity(c1, c2, metric=DataStructs.TanimotoSimilarity)
                    dice_sim = DataStructs.FingerprintSimilarity(c1, c2, metric=DataStructs.DiceSimilarity)
                    if dice_sim >= 0.5:
                        records[i]["similar_compounds"] = records[i].get("similar_compounds", []) + [ {"compound" : j, "morgan_fingerprint_2_dice":dice_sim, "morgan_fingerprint_2_tanimoto":tan_sim } ]

                        print(json.dumps( {
                            "compound_1" : i,
                            "compound_2" : j,
                            "morgan_fingerprint_2_dice":dice_sim, 
                            "morgan_fingerprint_2_tanimoto":tan_sim
                        } ))

    #for v in records.values():
    #    print(json.dumps(v))
