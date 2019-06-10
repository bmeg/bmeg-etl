import pandas

import bmeg.ioutils
from bmeg import (Sample, Aliquot, Case, Project, Program,
                  Aliquot_Sample_Sample,
                  Sample_Case_Case,
                  Case_Projects_Project,
                  Sample_Projects_Project,
                  Aliquot_Projects_Project,
                  Case_Phenotypes_Phenotype,
                  Sample_Phenotypes_Phenotype,
                  Project_Programs_Program)
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from bmeg.emitter import JSONEmitter


def transform(cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              phenotype_lookup_path="source/ccle/cellline_phenotype_lookup.tsv",
              cellline_meta_path="source/gdsc/Cell_Lines_Details.xlsx",
              drugs_meta_path="source/gdsc/Screened_Compounds.xlsx",
              emitter_prefix="gdsc",
              emitter_directory="gdsc"):

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
    projects = bmeg.ioutils.read_lookup(project_lookup_path)
    phenotypes = bmeg.ioutils.read_lookup(phenotype_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    prog = Program(submitter_id=Program.make_gid("GDSC"),
                   program_id="GDSC")
    emitter.emit_vertex(prog)

    drugs_df = pandas.read_excel(drugs_meta_path, sheet_name=0, header=0)
    cells_df = pandas.read_excel(cellline_meta_path, sheet_name=0, header=0)
    drugs = drugs_df["DRUG_NAME"].tolist()

    emitted_celllines = {}
    emitted_projects = {}
    for i, row in cells_df.iterrows():
        emit_cellline = False
        cosmic = row.get("COSMIC identifier")
        if pandas.isnull(cosmic):
            cosmic = None
        else:
            cosmic = str(int(cosmic))
        if cosmic in celllines:
            cellline_id = celllines[cosmic]
        elif str(row["Sample Name"]) in celllines:
            cellline_id = celllines[str(row["Sample Name"])]
        elif str(row["Sample Name"]).replace("-", "").replace("/", "").upper() in celllines:
            cellline_id = celllines[str(row["Sample Name"]).replace("-", "").replace("/", "").upper()]
        else:
            cellline_id = str(row["Sample Name"])
            emit_cellline = True

        if cellline_id == "TOTAL:":
            continue

        if cellline_id in emitted_celllines:
            continue

        project_id = "GDSC_%s" % (projects.get(cellline_id, "Unknown"))
        proj = Project(submitter_id=Project.make_gid(project_id), project_id=project_id)
        if proj.gid() not in emitted_projects:
            emitter.emit_vertex(proj)
            emitter.emit_edge(
                Project_Programs_Program(
                    from_gid=proj.gid(),
                    to_gid=prog.gid(),
                ),
                emit_backref=True
            )
            emitted_projects[proj.gid()] = None

        # cellline cases belong to multiple projects so we leave project_id blank...
        c = Case(submitter_id=Case.make_gid(cellline_id),
                 case_id=cellline_id,
                 project_id='')
        if emit_cellline:
            emitter.emit_vertex(c)
            # case <-> project edges
            emitter.emit_edge(
                Case_Projects_Project(
                    from_gid=c.gid(),
                    to_gid=proj.gid()
                ),
                emit_backref=True
            )

        sample_id = "GDSC:%s" % (cellline_id)
        s = Sample(submitter_id=Sample.make_gid(sample_id),
                   sample_id=sample_id,
                   project_id=proj.gid())
        emitter.emit_vertex(s)
        # sample <-> case edges
        emitter.emit_edge(
            Sample_Case_Case(
                from_gid=s.gid(),
                to_gid=c.gid()
            ),
            emit_backref=True
        )
        # sample <-> project edges
        emitter.emit_edge(
            Sample_Projects_Project(
                from_gid=s.gid(),
                to_gid=proj.gid()
            ),
            emit_backref=True
        )

        phenotype_name = phenotypes.get(cellline_id, None)
        if phenotype_name:
            pheno = phenotype_factory(phenotype_name)
            emitter.emit_vertex(pheno)
            # case <-> phenotype edges
            emitter.emit_edge(
                Case_Phenotypes_Phenotype(
                    from_gid=c.gid(),
                    to_gid=pheno.gid()
                ),
                emit_backref=True
            )
            # sample <-> phenotype edges
            emitter.emit_edge(
                Sample_Phenotypes_Phenotype(
                    from_gid=s.gid(),
                    to_gid=pheno.gid()
                ),
                emit_backref=True
            )

        experiement_type = "DrugResponse"
        for drug in drugs:
            aliquot_id = "GDSC:%s:%s:%s" % (cellline_id, experiement_type, drug)
            a = Aliquot(submitter_id=Aliquot.make_gid(aliquot_id),
                        aliquot_id=aliquot_id,
                        project_id=proj.gid())
            emitter.emit_vertex(a)
            # aliquot <-> sample edges
            emitter.emit_edge(
                Aliquot_Sample_Sample(
                    from_gid=a.gid(),
                    to_gid=s.gid()
                ),
                emit_backref=True
            )
            # aliquot <-> project edges
            emitter.emit_edge(
                Aliquot_Projects_Project(
                    from_gid=a.gid(),
                    to_gid=proj.gid()
                ),
                emit_backref=True
            )

        emitted_celllines[cellline_id] = None

    emitter.close()


if __name__ == "__main__":
    transform()
