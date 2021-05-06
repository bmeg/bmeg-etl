#!/usr/bin/env python

import argparse
import json
import sqlite3

from rdkit import Chem
from rdkit.Chem import AllChem


def genObjects(conn, ids):
    c = conn.cursor()

    for i in ids:
        c.execute("select * from MOLECULE_DICTIONARY where CHEMBL_ID == ?", [i])
        names = list(i[0] for i in c.description)
        entry = {}
        for row in c:
            entry = dict(zip(names,row))
        if "molregno" in entry:
            synonyms = []
            for row in c.execute("select distinct(SYNONYMS) from molecule_synonyms where molregno=?", [entry["molregno"]]):
                synonyms.append(row[0])
            entry['synonyms'] = synonyms
            smiles = None
            inchi = None
            inchiKey = None
            for row in c.execute("select STANDARD_INCHI, STANDARD_INCHI_KEY, CANONICAL_SMILES from COMPOUND_STRUCTURES where MOLREGNO=?", [entry["molregno"]]):
                inchi=row[0]
                inchiKey = row[1]
                smiles = row[2]
            if smiles is not None:
                entry["smiles"] = smiles
                m = Chem.MolFromSmiles(smiles)
                fp = AllChem.GetMorganFingerprintAsBitVect(m, radius=2)
                fingerprint = list(fp)
                entry["morgan_fingerprint_2"] = fingerprint

            if inchi is not None:
                entry["inchi"] = inchi
                entry["inchi_key"] = inchiKey
            entry["id"] = "Compound:%s" % (i)

            print(json.dumps(entry))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", type=int, default=0)
    parser.add_argument("db")
    parser.add_argument("input")

    args = parser.parse_args()
    conn = sqlite3.connect(args.db)
    ids = []
    with open(args.input) as handle:
        for line in handle:
            tmp = line.rstrip().split("\t")
            ids.append(tmp[args.c])
    genObjects(conn, ids)
