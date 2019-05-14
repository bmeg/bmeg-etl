import pandas

import bmeg.ioutils
from bmeg.vertex import Sample, Aliquot, Case, Project, Program
from bmeg.edge import AliquotFor, SampleFor, InProject, InProgram, PhenotypeOf
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from bmeg.emitter import JSONEmitter


def transform(cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              phenotype_lookup_path="source/ccle/cellline_phenotype_lookup.tsv",
              metadrugPath='source/ctrp/v20.meta.per_compound.txt',
              metacelllinePath="source/ctrp/v20.meta.per_cell_line.txt",
              responsePath="source/ctrp/v20.data.curves_post_qc.txt",
              metaexperimentPath="source/ctrp/v20.meta.per_experiment.txt",
              curvePath="source/ctrp/v20.data.per_cpd_post_qc.txt",
              emitter_prefix='ctrp',
              emitter_directory='ctrp'):

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
    projects = bmeg.ioutils.read_lookup(project_lookup_path)
    phenotypes = bmeg.ioutils.read_lookup(phenotype_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    prog = Program(program_id="CTRP")
    emitter.emit_vertex(prog)

    compound_df = pandas.read_table(metadrugPath)
    compound_df = compound_df.set_index("master_cpd_id")

    ccl_df = pandas.read_table(metacelllinePath)
    ccl_df = ccl_df.set_index("master_ccl_id")

    response_df = pandas.read_table(responsePath)

    metaexperiment_df = pandas.read_table(metaexperimentPath)
    metaexperiment_df = metaexperiment_df.drop_duplicates("experiment_id")
    metaexperiment_df = metaexperiment_df.set_index("experiment_id")

    curve_df = pandas.read_table(curvePath)
    curve_df = curve_df.set_index(["experiment_id", "master_cpd_id"])

    raw_ids = {}
    drugs = {}
    for i, row in response_df.iterrows():
        cpd_id = row['master_cpd_id']
        cpd_name = compound_df.loc[cpd_id]['cpd_name']
        exp_id = row['experiment_id']
        ccl_id = metaexperiment_df.loc[exp_id]['master_ccl_id']
        ccl_name = ccl_df.loc[ccl_id]['ccl_name']

        # both ids map to ACH-000991
        if ccl_name in ["697", "SNU81"]:
            ccl_name = "SNU81"

        if ccl_name not in raw_ids:
            raw_ids[ccl_name] = None
        if ccl_name in drugs:
            drugs[ccl_name].append(cpd_name)
        else:
            drugs[ccl_name] = [cpd_name]

    emitted_celllines = {}
    emitted_projects = {}
    for i in raw_ids:
        emit_cellline = False
        if i in celllines:
            cellline_id = celllines[i]
        elif i.split("_")[0] in celllines:
            cellline_id = celllines[i.split("_")[0]]
        else:
            emit_cellline = True
            cellline_id = i

        if cellline_id in emitted_celllines:
            continue

        project_id = "CTRP_%s" % (projects.get(cellline_id, "Unknown"))
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

        sample_id = "CTRP:%s" % (cellline_id)
        s = Sample(sample_id=sample_id)
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
        for drug in drugs.get(i, []):
            aliquot_id = "CTRP:%s:%s:%s" % (cellline_id, experiement_type, drug)
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


if __name__ == '__main__':  # pragma: no cover
    transform()
