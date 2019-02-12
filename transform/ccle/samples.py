import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample, Aliquot, Case, Project, Program
from bmeg.edge import AliquotFor, BiosampleFor, InProject, InProgram, PhenotypeOf
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from pydash import is_blank


def transform(path="source/ccle/DepMap-2018q4-celllines.csv",
              prefix="ccle"):
    emitter = JSONEmitter(directory=prefix)
    reader = bmeg.ioutils.read_csv(path)

    # Load sample metadata.
    case_gids = []
    project_gids = []
    phenotype_gids = []
    prog = Program(program_id="DepMap")
    emitter.emit_vertex(prog)
    # DepMap_ID,CCLE_Name,Aliases,COSMIC_ID,Sanger ID,Primary Disease,Subtype Disease,Gender,Source
    # ACH-000001,NIHOVCAR3_OVARY,NIH:OVCAR-3;OVCAR3,905933,2201,Ovarian Cancer,"Adenocarcinoma, high grade serous",Female,ATCC
    for row in reader:
        sample_id = row["DepMap_ID"]
        b = Biosample(biosample_id=sample_id,
                      ccle_attributes=row)
        emitter.emit_vertex(b)

        a = Aliquot(aliquot_id=sample_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid()
        )

        i = Case(case_id=sample_id,
                 ccle_attributes={'gender': row.get('Gender', None)})
        if i.gid() not in case_gids:
            emitter.emit_vertex(i)
            case_gids.append(i.gid())
        emitter.emit_edge(
            BiosampleFor(),
            b.gid(),
            i.gid()
        )

        project_id = "DepMap_{}".format("_".join(row["Primary Disease"].split()))
        # create project
        proj = Project(project_id=project_id)
        if proj.gid() not in project_gids:
            emitter.emit_vertex(proj)
            emitter.emit_edge(
                InProgram(),
                proj.gid(),
                prog.gid()
            )
            project_gids.append(proj.gid())
        emitter.emit_edge(
            InProject(),
            i.gid(),
            proj.gid()
        )

        phenotype_name = row.get('Subtype Disease', None)
        if not phenotype_name or is_blank(phenotype_name):
            phenotype_name = row.get('Primary Disease')
        phenotype = phenotype_factory(name=phenotype_name)
        if phenotype.gid() not in phenotype_gids:
            emitter.emit_vertex(phenotype)
            phenotype_gids.append(phenotype.gid())
        emitter.emit_edge(
            PhenotypeOf(),
            a.gid(),
            phenotype.gid(),
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
