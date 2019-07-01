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

    prog = Program(id=Program.make_gid("CTRP"),
                   program_id="CTRP")
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
    emitted_phenotypes = {}
    for i in raw_ids:
        if i in celllines:
            cellline_id = celllines[i]
        elif i.split("_")[0] in celllines:
            cellline_id = celllines[i.split("_")[0]]
        else:
            cellline_id = i

        if cellline_id in emitted_celllines:
            continue

        project_id = "CTRP_%s" % (projects.get(cellline_id, "Unknown"))
        proj = Project(id=Project.make_gid(project_id),
                       project_id=project_id)
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

        c = Case(id=Case.make_gid(cellline_id),
                 submitter_id=i,
                 case_id=cellline_id,
                 project_id=proj.gid())
        emitter.emit_vertex(c)
        # case <-> project edges
        emitter.emit_edge(
            Case_Projects_Project(
                from_gid=c.gid(),
                to_gid=proj.gid()
            ),
            emit_backref=True
        )

        sample_id = "CTRP:%s" % (cellline_id)
        s = Sample(id=Sample.make_gid(sample_id),
                   submitter_id=i,
                   sample_id=cellline_id,
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
            # sample <-> phenotype edges
            emitter.emit_edge(
                Sample_Phenotypes_Phenotype(
                    from_gid=s.gid(),
                    to_gid=pheno.gid()
                ),
                emit_backref=True
            )

        experiement_type = "DrugResponse"
        for drug in drugs.get(i, []):
            aliquot_id = "CTRP:%s:%s:%s" % (cellline_id, experiement_type, drug)
            a = Aliquot(id=Aliquot.make_gid(aliquot_id),
                        submitter_id=i,
                        aliquot_id=cellline_id,
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


if __name__ == '__main__':  # pragma: no cover
    transform()
