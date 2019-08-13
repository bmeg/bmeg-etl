import logging
import bmeg.enrichers.gene_enricher as gene_enricher

from bmeg import (Pathway, PathwayComponent, Interaction, Gene, Publication, Project,
                  Pathway_PathwayComponents_PathwayComponent, Pathway_Interactions_Interaction, Pathway_Genes_Gene, Pathway_SubPathway_Pathway,
                  PathwayComponent_Interactions_Interaction, PathwayComponent_Genes_Gene,
                  Interaction_Genes_Gene, Interaction_Publications_Publication,
                  Gene_PathwayComponents_PathwayComponent, Gene_Interactions_Interaction)
from bmeg.emitter import JSONEmitter
from bmeg.ioutils import read_tsv


def get_participant(p):
    v = None
    try:
        gene = gene_enricher.get_gene(p)
        ens_id = gene.get("ensembl_gene_id", None)
        if ens_id is None:
            raise ValueError("No ensembl id found for entrez_id: {}".format(p))
        v = Gene(id=Gene.make_gid(ens_id))
    except Exception as e:
        logging.error(e)
    return v


def transform(sif_file="source/pathway_commons/pc11.details.sif",
              pathways_file="source/pathway_commons/pathways.txt",
              emitter_prefix=None,
              emitter_directory="pathway_commons"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    sif = read_tsv(sif_file,
                   fieldnames=["PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B",
                               "MEDIATOR", "PUBMED", "PMC", "COMMENTS", "PATHWAY",
                               "PATHWAY_URI", "RESOURCE", "SOURCE_LOC", "TARGET_LOC"])

    components = {}
    pathways = {}
    path_comp = {}
    path_int = {}
    gene_comp = {}
    gene_path = {}
    for line in sif:
        # ignore compounds production for now
        if line["PARTICIPANT_A"].startswith("CHEBI") or line["PARTICIPANT_B"].startswith("CHEBI"):
            continue

        i = Interaction(
            id=Interaction.make_gid(line["PARTICIPANT_A"], line["INTERACTION_TYPE"], line["PARTICIPANT_B"]),
            Project=Project.make_gid("Reference"),
            source=line["RESOURCE"],
            type=line["INTERACTION_TYPE"]
        )
        pa = get_participant(line["PARTICIPANT_A"])
        pb = get_participant(line["PARTICIPANT_B"])

        # emit mapped interactions
        if pa is None or pb is None:
            continue
        emitter.emit_vertex(i)
        emitter.emit_edge(
            Gene_Interactions_Interaction(
                from_gid=pa.gid(),
                to_gid=i.gid()
            ),
            emit_backref=False
        )
        emitter.emit_edge(
            Interaction_Genes_Gene(
                from_gid=i.gid(),
                to_gid=pb.gid()
            ),
            emit_backref=False
        )

        # emit publication edges
        for p in line["PUBMED"].split(";"):
            if p != "":
                emitter.emit_edge(
                    Interaction_Publications_Publication(
                        from_gid=i.gid(),
                        to_gid=Publication.make_gid("ncbi.nlm.nih.gov/pubmed/{}".format(p))
                    ),
                    emit_backref=True
                )

        # emit pathway components and edges
        for m in line["MEDIATOR"].split(";"):
            comp = PathwayComponent(
                id=PathwayComponent.make_gid(m),
                Project=Project.make_gid("Reference")
            )
            if m not in components:
                emitter.emit_vertex(comp)
                components[m] = True
            emitter.emit_edge(
                PathwayComponent_Interactions_Interaction(
                    from_gid=comp.gid(),
                    to_gid=i.gid()
                ),
                emit_backref=True
            )
            if (pa.gid(), comp.gid()) not in gene_comp:
                emitter.emit_edge(
                    Gene_PathwayComponents_PathwayComponent(
                        from_gid=pa.gid(),
                        to_gid=comp.gid()
                    ),
                    emit_backref=False
                )
                gene_comp[(pa.gid(), comp.gid())] = True
            if (comp.gid(), pb.gid()) not in gene_comp:
                emitter.emit_edge(
                    PathwayComponent_Genes_Gene(
                        from_gid=comp.gid(),
                        to_gid=pb.gid()
                    ),
                    emit_backref=False
                )
                gene_comp[(comp.gid(), pb.gid())] = True

            # emit pathway edges
            for p in line["PATHWAY_URI"].split(";"):
                if p != "":
                    pathway = Pathway(
                        id=Pathway.make_gid(p),
                        Project=Project.make_gid("Reference"),
                        name=line["PATHWAY"]
                    )
                    if pathway.gid() not in pathways:
                        emitter.emit_vertex(pathway)
                        pathways[pathway.gid()] = True
                    if comp.gid() not in path_comp:
                        emitter.emit_edge(
                            Pathway_PathwayComponents_PathwayComponent(
                                from_gid=pathway.gid(),
                                to_gid=comp.gid(),
                            ),
                            emit_backref=True
                        )
                        path_comp[comp.gid()] = True
                    if i.gid() not in path_int:
                        emitter.emit_edge(
                            Pathway_Interactions_Interaction(
                                from_gid=pathway.gid(),
                                to_gid=i.gid(),
                            ),
                            emit_backref=True
                        )
                        path_int[i.gid()] = True
                    # Pathway -> Gene edges are bidirectional
                    if pa.gid() not in gene_path:
                        emitter.emit_edge(
                            Pathway_Genes_Gene(
                                from_gid=pathway.gid(),
                                to_gid=pa.gid()
                            ),
                            emit_backref=True
                        )
                        gene_path[pa.gid()] = True
                    if pb.gid() not in gene_path:
                        emitter.emit_edge(
                            Pathway_Genes_Gene(
                                from_gid=pathway.gid(),
                                to_gid=pb.gid()
                            ),
                            emit_backref=True
                        )
                        gene_path[pb.gid()] = True

    # emit pathway hierarchy
    # PATHWAY_URI	DISPLAY_NAME	DIRECT_SUB_PATHWAY_URIS	ALL_SUB_PATHWAY_URIS
    sub_paths = {}
    pathways = read_tsv(pathways_file)
    for line in pathways:
        pathway_gid = Pathway.make_gid(line["PATHWAY_URI"])
        for sub in line["ALL_SUB_PATHWAY_URIS"].split(";"):
            sub_gid = Pathway.make_gid(line[sub])
            if sub_gid not in sub_paths:
                emitter.emit_edge(
                    Pathway_SubPathway_Pathway(
                        from_gid=pathway_gid,
                        to_gid=sub_gid,
                    ),
                    emit_backref=True
                )
                sub_paths[sub_gid] = True

    emitter.close()


if __name__ == "__main__":
    transform()
