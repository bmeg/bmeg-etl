from bmeg.models.ensembl_lookup import hugo_ensembl

import mygene


def query_mygeneinfo_for_ensembl_gene(id):
    """
    Query http://mygene.info/v3/api for ensembl gene id
    """
    if not isinstance(id, str):
        raise TypeError("expected str")

    mg = mygene.MyGeneInfo()
    res = mg.querymany(
        id,
        scopes="accession,alias,ensembl.gene,ensembl.transcript,ensembl.protein,\
        entrezgene,name,refseq,symbol,unigene,uniprot",
        fields="ensembl.gene",
        species="human"
    )

    # TODO
    # is it ok to always take the highest scoring hit?
    top_hit = res[0]
    if 'notfound' in top_hit or "ensembl" not in top_hit:
        raise RuntimeError("query for %s returned no hits" % id)

    # TODO
    # should we always pick the first item if 'ensembl' is a list?
    if isinstance(top_hit["ensembl"], list):
        top_hit["ensembl"] = top_hit["ensembl"][0]

    return top_hit["ensembl"]["gene"]


def ensembl_gene_lookup(id):
    ensembl_id = hugo_ensembl(id)
    if ensembl_id == "":
        raise ValueError("ensembl id for %s not found" % id)
    return ensembl_id


def ensembl_transcript_lookup(id):
    raise RuntimeError("not implemented")


def ensembl_protein_lookup(id):
    raise RuntimeError("not implemented")


def disease_ontology_lookup(id):
    raise RuntimeError("not implemented")
