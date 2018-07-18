import ensembl_lookup


def ensembl_gene_lookup(id):
    ensembl_id = ensembl_lookup.hugo_ensembl(id)
    if ensembl_id == "":
        raise ValueError("ensembl id for %s not found" % id)
    return ensembl_id


def ensembl_transcript_lookup(id):
    raise RuntimeError("not implemented")


def ensembl_protein_lookup(id):
    raise RuntimeError("not implemented")


def disease_ontology_lookup(id):
    raise RuntimeError("not implemented")
