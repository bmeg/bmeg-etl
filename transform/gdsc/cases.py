import pandas

import bmeg.ioutils
from bmeg.vertex import Sample, Aliquot, Case, Project, Program
from bmeg.edge import AliquotFor, SampleFor, InProject, InProgram, PhenotypeOf
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

    prog = Program(program_id="GDSC")
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

        if cellline_id in emitted_celllines:
            continue

        project_id = "GDSC_%s" % (projects.get(cellline_id, "Unknown"))
        p = Project(project_id=project_id)
        if p.gid() not in emitted_projects:
            emitter.emit_vertex(p)
            emitter.emit_edge(
                InProgram(),
                p.gid(),
                prog.gid(),
            )
            emitted_projects[p.gid()] = None

        c = Case(case_id=cellline_id)
        if emit_cellline:
            emitter.emit_vertex(c)
        emitter.emit_edge(
            InProject(),
            c.gid(),
            p.gid(),
        )

        s = Sample(sample_id=cellline_id)
        emitter.emit_vertex(s)
        emitter.emit_edge(
            SampleFor(),
            s.gid(),
            c.gid(),
        )
        emitter.emit_edge(
            InProject(),
            s.gid(),
            p.gid(),
        )

        phenotype_name = phenotypes.get(cellline_id, None)
        if phenotype_name:
            pheno = phenotype_factory(phenotype_name)
            emitter.emit_vertex(pheno)
            emitter.emit_edge(
                PhenotypeOf(),
                c.gid(),
                pheno.gid()
            )
            emitter.emit_edge(
                PhenotypeOf(),
                s.gid(),
                pheno.gid()
            )

        experiement_type = "DrugResponse"
        for drug in drugs:
            aliquot_id = "%s:%s:%s" % (cellline_id, experiement_type, drug)
            a = Aliquot(aliquot_id=aliquot_id)
            emitter.emit_vertex(a)
            emitter.emit_edge(
                AliquotFor(),
                a.gid(),
                s.gid(),
            )
            emitter.emit_edge(
                InProject(),
                a.gid(),
                p.gid(),
            )

        emitted_celllines[cellline_id] = None

    emitter.close()


if __name__ == "__main__":
    transform()
