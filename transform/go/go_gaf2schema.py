#!/usr/bin/env python

# parses GO annotation file
# http://www.geneontology.org/doc/GO.references

import gzip

from bmeg import Gene, GeneOntologyTerm, GeneOntologyTerm_Genes_Gene
from bmeg.emitter import JSONEmitter


UNIPROT_COL = 1
SYMBOL_COL = 2
GOID_COL = 4
REF_COL = 5
EVIDENCE_COL = 6
NAME_COL = 9


def transform(gaf_file="source/go/goa_human.gaf.gz",
              id_map_file="source/go/HUMAN_9606_idmapping.dat.gz",
              emitter_prefix=None,
              emitter_directory="go"):

    uniprot_2_ensembl = {}

    with gzip.GzipFile(id_map_file) as handle:
        for line in handle:
            row = line.decode().rstrip().split("\t")
            if row[1] == "Ensembl":
                uniprot_2_ensembl[row[0]] = row[2]

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    with gzip.GzipFile(gaf_file) as handle:
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

                    gene_gid = Gene.make_gid(ensembl_id)
                    go_gid = GeneOntologyTerm.make_gid(go_id)

                    emitter.emit_edge(
                        GeneOntologyTerm_Genes_Gene(
                            to_gid=gene_gid,
                            from_gid=go_gid,
                            data={
                                "evidence": evidence,
                                "references": references,
                                "title": title
                            }
                        ),
                        emit_backref=True
                    )

    emitter.close()


if __name__ == "__main__":
    transform()
