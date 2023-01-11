#!/usr/bin/env python

import sys
import json
from rdkit import Chem
from rdkit.Chem import AllChem


for line in sys.stdin:
    row = json.loads(line)
    if "canonical_smiles" in row:
        smiles = row["canonical_smiles"]
        m = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(m, radius=2)
        fingerprint = list(fp)
        row["morgan_fingerprint_2"] = fingerprint
    print(json.dumps(row))
