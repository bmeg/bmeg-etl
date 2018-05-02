#!/usr/bin/env python

# parses file ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz

import sys
import gzip
import itertools
from bmeg.protein_pb2 import Protein

def parse_id(line):
    row = line.split(" ")
    return "uniprot_id", row[0]

def parse_rx(line):
    row = line.split("; ")
    out = []
    for i in row:
        if i.startswith("PubMed="):
            out.append(i.replace("PubMed=", ""))
    return "pubmed", out

def parse_ac(line):
    return "accession", line.rstrip("\n;").split("; ")

def parse_dr(line):
    row = line.split("; ")
    if row[0] == "PDB":
        return "pdb", row[1]
    if row[0] == "Ensembl":
        return "ensembl", [
            ("transcript", row[1]),
            ("protein", row[2]),
            ("gene", row[3].split(".")[0]),
        ]
    return None, None

def parse_de(line):
    if line.startswith("RecName: Full="):
        return "name", line.split("Full=")[1]
    return None, None

parse_map = {
    "ID" : parse_id,
    "RX" : parse_rx,
    "AC" : parse_ac,
    "DR" : parse_dr,
    "DE" : parse_de
}

def message_prep(doc):
    out = Protein()
    out.uniprot_id = doc["uniprot_id"][0]
    out.name = doc["name"][0]
    if "pubmed" in doc:
        for i in list(itertools.chain.from_iterable(doc["pubmed"])):
            out.pubmed.append(i)
    for i in list(itertools.chain.from_iterable(doc["accession"])):
        out.accession.append(i)
    if 'pdb' in doc:
        for i in list(itertools.chain(doc["pdb"])):
            out.pdb.append(i)
    if 'ensembl' in doc:
        for ent in doc['ensembl']:
            for k, v in ent:
                if k == "transcript":
                    out.ensembl_transcript = v
                if k == "protein":
                    out.ensembl_protein = v
                if k == "gene":
                    out.ensembl_gene = v
    return out

def sprot_parse(path):
    cur = {}
    with gzip.GzipFile(path) as handle:
        for line in handle:
            key = line[:2]
            value = line.rstrip()[5:]
            
            if key == "//":
                print message_prep(cur)
                cur = {}
            if key in parse_map:
                k,v = parse_map[key](value)
                if v is not None:
                    cur[k] = cur.get(k, []) + [v]
            


if __name__ == "__main__":
    
    sprot_parse(sys.argv[1])