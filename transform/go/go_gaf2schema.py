#!/usr/bin/env python

# parses GO annotation file
# http://www.geneontology.org/doc/GO.references


import sys
import json
import gzip

from bmeg.vertex import Gene, GeneOntologyTerm
from bmeg.edge import GeneOntologyAnnotation
from bmeg.emitter import JSONEmitter


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
