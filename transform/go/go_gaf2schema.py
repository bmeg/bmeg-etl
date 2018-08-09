#!/usr/bin/env python

# parses GO annotation file
# http://www.geneontology.org/doc/GO.references

# curl -L -O http://geneontology.org/gene-associations/goa_human.gaf.gz
# curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz

import sys
import json
import gzip

from bmeg.vertex import Gene, GeneOntologyTerm
from bmeg.edge import GeneOntologyAnnotation
from bmeg.emitter import JSONEmitter

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

	with gzip.GzipFile(sys.argv[2]) as handle:
		for line in handle:
			row = line.decode().rstrip().split("\t")
			if row[1] == "Ensembl":
				uniprot_2_ensembl[row[0]] = row[2]

	emitter = JSONEmitter(sys.argv[3])
	with gzip.GzipFile(sys.argv[1]) as handle:
		for line in handle:
			line = line.decode()
			if not line.startswith("!"):
				row = line.rstrip().rsplit("\t")
				if row[UNIPROT_COL] in uniprot_2_ensembl:
					go_id = row[GOID_COL]
					#for r in row[SYMBOL_COL].split("|"):
					#	if re.search(r'^[A-Z0-9]+$', r):
					#		goa.genes.append( gene_ns(r) )
					ensembl_id = uniprot_2_ensembl[row[UNIPROT_COL]]
					evidence = row[EVIDENCE_COL]
					references = []
					if row[REF_COL].startswith("GO_REF:"):
						references.append(row[REF_COL])
					if row[REF_COL].startswith("PMID:"):
						references.append(row[REF_COL])

					if len(row[NAME_COL]):
						title = row[NAME_COL]

					gene_gid = Gene.make_gid(gene_id=ensembl_id)
					go_gid = GeneOntologyTerm.make_gid(go_id=go_id)
					emitter.emit_edge(GeneOntologyAnnotation(
							evidence=evidence,
							references=references, title=title),
						to_gid=gene_gid,
						from_gid=go_gid
					)
	emitter.close()
