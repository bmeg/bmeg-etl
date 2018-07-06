#!/usr/bin/env python

# parses GO annotation file
# http://www.geneontology.org/doc/GO.references

# curl -L -O http://geneontology.org/gene-associations/goa_human.gaf.gz
# curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz

import sys
import csv
import re
import json
import gzip
from bmeg import phenotype_pb2
from google.protobuf import json_format

def uniprot_ns(n):
	return "uniprot:" + n

def ensembl_ns(n):
	return "ensembl:" + n

def go_ns(n):
	return "go:" + n

def goref_ns(n):
	return "goref:" + n

def gene_ns(n):
	return "gene:" + n

def pubmed_ns(n):
	return "pubmed:" + n

UNIPROT_COL = 1
SYMBOL_COL = 2
GOID_COL = 4
REF_COL = 5
EVIDENCE_COL = 6
NAME_COL = 9

def message_to_json(handle, message):
	msg = json_format.MessageToDict(message, preserving_proto_field_name=True)
	handle.write(json.dumps(msg))
	handle.write("\n")

if __name__ == "__main__":

	uniprot_2_ensembl = {}

	output = open(sys.argv[3], "w")

	with gzip.GzipFile(sys.argv[2]) as handle:
		for line in handle:
			row = line.rstrip().split("\t")
			if row[1] == "Ensembl":
				uniprot_2_ensembl[row[0]] = row[2]


	with gzip.GzipFile(sys.argv[1]) as handle:
		for line in handle:
			if not line.startswith("!"):
				row = line.rstrip().rsplit("\t")
				if row[UNIPROT_COL] in uniprot_2_ensembl:
					go_id = row[GOID_COL]
					goa = phenotype_pb2.GeneOntologyAnnotation()
					#for r in row[SYMBOL_COL].split("|"):
					#	if re.search(r'^[A-Z0-9]+$', r):
					#		goa.genes.append( gene_ns(r) )
					goa.genes.append(uniprot_2_ensembl[row[UNIPROT_COL]])
					goa.functions.append(go_id)
					goa.evidence.append(row[EVIDENCE_COL])

					if row[REF_COL].startswith("GO_REF:"):
						goa.references.append(row[REF_COL])
					if row[REF_COL].startswith("PMID:"):
						goa.references.append(row[REF_COL])

					if len(row[NAME_COL]):
						goa.title = row[NAME_COL]

					message_to_json(output, goa)
	output.close()
