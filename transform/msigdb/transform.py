import logging
import xml.etree.ElementTree as ET
import bmeg.enrichers.gene_enricher as gene_enricher

from bmeg import (GeneSet, Gene, Publication, Project,
                  GeneSet_Genes_Gene, GeneSet_Publications_Publication)
from bmeg.emitter import JSONEmitter


def transform(input_path="source/msigdb/msigdb_v6.2.xml",
              emitter_prefix=None,
              emitter_directory="msigdb"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    root = ET.parse(input_path).getroot()
    # all props
    # ['standard_name', 'systematic_name', 'historical_names', 'organism', 'pmid',
    #  'authors', 'geoid', 'exact_source', 'geneset_listing_url', 'external_details_url',
    #  'chip', 'category_code', 'sub_category_code', 'contributor', 'contributor_org',
    #  'description_brief', 'description_full', 'tags', 'members', 'members_symbolized',
    #  'members_ezid', 'members_mapping', 'founder_names', 'refinement_datasets',
    #  'validation_datasets']
    keep_props = [
        'standard_name', 'systematic_name', 'historical_names', 'geoid', 'exact_source',
        'geneset_listing_url', 'external_details_url', 'chip', 'category_code', 'sub_category_code',
        'contributor', 'contributor_org', 'description_brief', 'description_full'
    ]
    for gs in root.findall("GENESET"):
        props = {}
        for k, v in gs.items():
            if k.lower() in keep_props:
                props[k.lower()] = v

        # skip pathways since those are defined elsewhere
        if props["category_code"] == "C2" and props["sub_category_code"] == "CP":
            continue

        # skip data with license restrictions
        # http://software.broadinstitute.org/gsea/msigdb_license_terms.jsp
        if props["standard_name"].startswith("BIOCARTA_") or \
           props["standard_name"].startswith("ST_") or \
           props["standard_name"].startswith("KEGG_"):
            continue

        # skip gene ontologies since those are defined elsewhere
        if props["category_code"] == "C5":
            continue

        props['id'] = GeneSet.make_gid(props["standard_name"])
        props['project_id'] = Project.make_gid("Reference")
        gene_set = GeneSet(**props)
        emitter.emit_vertex(gene_set)

        # create edge(s) to publications
        pmid = gs.get("PMID", "")
        if pmid != "":
            emitter.emit_edge(
                GeneSet_Publications_Publication(
                    from_gid=gene_set.gid(),
                    to_gid=Publication.make_gid("ncbi.nlm.nih.gov/pubmed/{}".format(pmid))
                ),
                emit_backref=True
            )

        # create edges to genes
        members = gs.get("MEMBERS_EZID", "")
        # members = gs.get("MEMBERS_SYMBOLIZED", "")
        if members == "":
            raise Exception("No members found for {}".format(gene_set.gid()))

        for g in members.split(","):
            try:
                gene = gene_enricher.get_gene(g)
                ens_id = gene.get("ensembl_gene_id", None)
                if ens_id is None:
                    raise ValueError("No ensembl id found for entrez_id: {}".format(g))
                emitter.emit_edge(
                    GeneSet_Genes_Gene(
                        from_gid=gene_set.gid(),
                        to_gid=Gene.make_gid(ens_id)
                    ),
                    emit_backref=True
                )
            except Exception as e:
                logging.error(e)

    emitter.close()


if __name__ == "__main__":
    transform()
