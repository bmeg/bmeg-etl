import json
import logging

from bmeg.ioutils import read_tsv
from bmeg.emitter import JSONEmitter
from bmeg import (Gene, G2PAssociation, Project, Publication, Compound,
                  G2PAssociation_Compounds_Compound, G2PAssociation_Genes_Gene,
                  G2PAssociation_Publications_Publication)
from bmeg.enrichers.drug_enricher import normalize
import bmeg.enrichers.gene_enricher as gene_enricher


def transform(interactions_file="source/dgidb/interactions.tsv",
              emitter_prefix=None,
              emitter_directory="dgidb"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    normalized_compounds = {}
    emitted_compounds = {}
    interactions = read_tsv(interactions_file)
    # gene_name gene_claim_name entrez_id interaction_claim_source interaction_types drug_claim_name drug_claim_primary_name drug_name drug_chembl_id PMIDs
    for line in interactions:
        source = line["interaction_claim_source"]
        # remove associations already brought in by VICC G2P
        if source in ["CGI", "CIViC", "OncoKB", "CKB"]:
            continue
        assoc_params = {
            "source": source,
            "source_document": json.dumps(line),
            "description": "{} {} {}".format(line["drug_chembl_id"], line["interaction_types"], line["gene_name"]),
            "evidence_label": None,
            "response_type": line["interaction_types"],
            "oncogenic": None,
            "source_url": None
        }
        assoc = G2PAssociation(
            id=G2PAssociation.make_gid(**assoc_params),
            project_id=Project.make_gid("Reference"),
            **assoc_params
        )
        gene_name = None
        if line["entrez_id"] is not None and line["entrez_id"] != "":
            gene_name = line["entrez_id"]
        elif line["gene_name"] is not None and line["gene_name"] != "":
            gene_name = line["gene_name"]
        else:
            gene_name = line["gene_claim_name"]
        try:
            gene = gene_enricher.get_gene(gene_name)
            ens_id = gene.get("ensembl_gene_id", None)
            if ens_id is None:
                raise ValueError("No ensembl id found for: {}".format(gene_name))
        except Exception as e:
            logging.error(e)
            continue
        chem_name = None
        if line["drug_chembl_id"] is not None and line["drug_chembl_id"] != "":
            chem_name = line["drug_chembl_id"]
        elif line["drug_name"] is not None and line["drug_name"] != "":
            chem_name = line["drug_name"]
        else:
            chem_name = line["drug_claim_name"]
        if chem_name in normalized_compounds:
            cinfo = normalized_compounds[chem_name]
        else:
            cinfo = normalize(name=chem_name)
            normalized_compounds[chem_name] = cinfo
        # excluding compounds we couldn't match since the interactions file contains things
        # that obviously are not referencing specific compounds
        if cinfo is None:
            continue
        compound = Compound(
            submitter_id=chem_name,
            project_id=Project.make_gid("Reference"),
            **cinfo
        )
        compound.id = Compound.make_gid(cinfo["id"])
        if compound.gid() not in emitted_compounds:
            emitter.emit_vertex(compound)
            emitted_compounds[compound.gid()] = True
        emitter.emit_vertex(assoc)
        emitter.emit_edge(
            G2PAssociation_Compounds_Compound(
                from_gid=assoc.gid(),
                to_gid=compound.gid()
            ),
            emit_backref=True
        )
        emitter.emit_edge(
            G2PAssociation_Genes_Gene(
                from_gid=assoc.gid(),
                to_gid=Gene.make_gid(ens_id)
            ),
            emit_backref=True
        )
        if line["PMIDs"] is None or line["PMIDs"] != "":
            continue
        pubs = line["PMIDs"].split(",")
        for p in pubs:
            emitter.emit_edge(
                G2PAssociation_Publications_Publication(
                    from_gid=assoc.gid(),
                    to_gid=Publication.make_gid("ncbi.nlm.nih.gov/pubmed/{}".format(p))
                ),
                emit_backref=True
            )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
