import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import (Case, Project, Case_Phenotypes_Phenotype)
from bmeg.enrichers.phenotype_enricher import phenotype_factory


def transform(path="source/ccle/DepMap-2019q1-celllines.csv_v2.csv",
              phenotype_lookup_path="source/ccle/cellline_phenotype_lookup.tsv",
              emitter_prefix="depmap",
              emitter_directory="ccle"):

    phenotypes = bmeg.ioutils.read_lookup(phenotype_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    reader = bmeg.ioutils.read_csv(path)

    # Load sample metadata.
    case_gids = {}
    emitted_phenotypes = {}
    # DepMap_ID,CCLE_Name,Aliases,COSMIC_ID,Sanger ID,Primary Disease,Subtype Disease,Gender,Source
    # ACH-000001,NIHOVCAR3_OVARY,NIH:OVCAR-3;OVCAR3,905933,2201,Ovarian Cancer,"Adenocarcinoma, high grade serous",Female,ATCC
    for row in reader:
        if "MERGED" in row["CCLE_Name"]:
            continue

        props = {}
        for k, v in row.items():
            props["_".join(k.split())] = v

        c = Case(
            submitter_id=Case.make_gid(row["DepMap_ID"]),
            case_id=row["DepMap_ID"],
            cellline_attributes=props,
            project_id=Project.make_gid('Shared')
        )
        if c.gid() not in case_gids:
            emitter.emit_vertex(c)
            case_gids[c.gid()] = None

        phenotype_name = phenotypes.get(row["DepMap_ID"], None)
        if phenotype_name:
            pheno = phenotype_factory(phenotype_name)
            if pheno.gid() not in emitted_phenotypes:
                emitter.emit_vertex(pheno)
                emitted_phenotypes[pheno.gid()] = None
            # case <-> phenotype edges
            emitter.emit_edge(
                Case_Phenotypes_Phenotype(
                    from_gid=c.gid(),
                    to_gid=pheno.gid()
                ),
                emit_backref=True
            )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
