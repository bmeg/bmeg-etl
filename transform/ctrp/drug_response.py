import pandas

import bmeg.ioutils
from bmeg import (Aliquot, DrugResponse, Project,
                  Compound_Projects_Project,
                  DrugResponse_Aliquot_Aliquot,
                  DrugResponse_Compounds_Compound)
from bmeg.enrichers.drug_enricher import compound_factory
from bmeg.emitter import JSONEmitter


def transform(cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              metadrugPath="source/ctrp/v20.meta.per_compound.txt",
              metacelllinePath="source/ctrp/v20.meta.per_cell_line.txt",
              responsePath="source/ctrp/v20.data.curves_post_qc.txt",
              metaexperimentPath="source/ctrp/v20.meta.per_experiment.txt",
              curvePath="source/ctrp/v20.data.per_cpd_post_qc.txt",
              emitter_prefix="drug_response",
              emitter_directory="ctrp"):

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
    projects = bmeg.ioutils.read_lookup(project_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

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

    emitted_compounds = {}
    project_compounds = {}
    for i, row in response_df.iterrows():
        exp_id = row["experiment_id"]
        cpd_id = row["master_cpd_id"]
        ccl_id = metaexperiment_df.loc[exp_id]["master_ccl_id"]

        # remap sample name if possible
        ccl_name = ccl_df.loc[ccl_id]["ccl_name"]
        if ccl_name in celllines:
            ccl_name = celllines[ccl_name]

        # Track drugs for project
        project_id = "CTRP_%s" % (projects.get(ccl_name, "Unknown"))
        proj = Project(submitter_id=Project.make_gid(project_id),
                       project_id=project_id)
        if proj.gid() not in project_compounds:
            project_compounds[proj.gid()] = {}

        cpd_name = compound_df.loc[cpd_id]["cpd_name"]
        auc = row["area_under_curve"]
        ec50 = row["apparent_ec50_umol"]

        curve_sub = curve_df.loc[(exp_id, cpd_id)]
        conc = curve_sub["cpd_conc_umol"]
        resp = curve_sub["cpd_avg_pv"]

        # create drug response vertex
        dr = DrugResponse(submitter_id=DrugResponse.make_gid("CTRP", ccl_name, cpd_name),
                          submitter_compound_id=cpd_name,
                          auc=auc,
                          ec50=ec50,
                          doses_um=list(conc),
                          activity_data_median=list(resp),
                          project_id=proj.gid())
        emitter.emit_vertex(dr)
        #  and edge to aliquot
        emitter.emit_edge(
            DrugResponse_Aliquot_Aliquot(
                from_gid=dr.gid(),
                to_gid=Aliquot.make_gid("CTRP:%s:DrugResponse:%s" % (dr.submitter_id, dr.submitter_compound_id))
            ),
            emit_backref=True
        )

        # create an edge to compound
        compound = compound_factory(name=cpd_name)
        if compound.gid() not in emitted_compounds:
            emitter.emit_vertex(compound)
            emitted_compounds[compound.gid()] = None
        emitter.emit_edge(
            DrugResponse_Compounds_Compound(
                from_gid=dr.gid(),
                to_gid=compound.gid()
            ),
            emit_backref=True
        )

        # create edge from compound to project
        if compound.gid() not in project_compounds[proj.gid()]:
            emitter.emit_edge(
                Compound_Projects_Project(
                    from_gid=compound.gid(),
                    to_gid=proj.gid()
                ),
                emit_backref=True
            )
            project_compounds[proj.gid()][compound.gid()] = True

    emitter.close()


if __name__ == "__main__":  # pragma: no cover
    transform()
