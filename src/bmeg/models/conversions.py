from bmeg.models.ensembl_lookup import hugo_ensembl

import mygene


def query_mygeneinfo_for_ensembl_gene(ids):
    """
    Query http://mygene.info/v3/api for ensembl gene id
    """
    mg = mygene.MyGeneInfo()
    res = mg.querymany(
        ids,
        scopes="accession,alias,ensembl.gene,ensembl.transcript,ensembl.protein,\
        entrezgene,name,refseq,symbol,unigene,uniprot",
        fields="ensembl.gene",
        species="human"
    )
    if 'notfound' in res[0]:
        raise RuntimeError("query for %s returned no hits" % id)

    # TODO
    # figure out how to handle the case where query returns "ensembl" as a list
    if isinstance(res[0]["ensembl"], list):
        res[0]["ensembl"] = res[0]["ensembl"][0]

    # TODO
    # is it ok to always take the highest scoring hit?
    if isinstance(ids, str):
        return res[0]["ensembl"]["gene"]
    else:
        return [x["ensembl"]["gene"] for x in res]


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
