import json
import logging

from bmeg.ioutils import read_tsv
from bmeg.emitter import JSONEmitter
from bmeg import (Gene, G2PAssociation, Project, Publication,
                  G2PAssociation_Compounds_Compound, G2PAssociation_Genes_Gene,
                  G2PAssociation_Publications_Publication)
from bmeg.enrichers.drug_enricher import compound_factory
import bmeg.enrichers.gene_enricher as gene_enricher


def transform(interactions_file="source/dgidb/interactions.tsv",
              emitter_prefix=None,
              emitter_directory="dgidb"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    interactions = read_tsv(interactions_file)
    # gene_name gene_claim_name entrez_id interaction_claim_source interaction_types drug_claim_name drug_claim_primary_name drug_name drug_chembl_id PMIDs
    for line in interactions:
        print(line)
        assoc_params = {
            "source": line["interaction_claim_source"],
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
        try:
            gene = gene_enricher.get_gene(line["entrez_id"])
            ens_id = gene.get("ensembl_gene_id", None)
            if ens_id is None:
                raise ValueError("No ensembl id found for entrez_id: {}".format(line["entrez_id"]))
        except Exception as e:
            logging.error(e)
            continue
        compound = compound_factory(name=line["drug_chembl_id"])
        emitter.emit_vertex(compound)
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
