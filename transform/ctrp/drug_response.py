import pandas

import bmeg.ioutils
from bmeg.vertex import DrugResponse, Aliquot, Project
from bmeg.edge import ResponseIn, ResponseTo, TestedIn
from bmeg.enrichers.drug_enricher import compound_factory
from bmeg.emitter import JSONEmitter


def transform(cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              metadrugPath='source/ctrp/v20.meta.per_compound.txt',
              metacelllinePath="source/ctrp/v20.meta.per_cell_line.txt",
              responsePath="source/ctrp/v20.data.curves_post_qc.txt",
              metaexperimentPath="source/ctrp/v20.meta.per_experiment.txt",
              curvePath="source/ctrp/v20.data.per_cpd_post_qc.txt",
              emitter_prefix='ctrp',
              emitter_directory='ctrp'):

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
    projects = bmeg.ioutils.read_lookup(project_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    compound_df = pandas.read_table(metadrugPath)
    compound_df = compound_df.set_index("master_cpd_id")
    for i, row in compound_df.iterrows():
        cpd_name = row['cpd_name']
        compound = compound_factory(name=cpd_name)
        emitter.emit_vertex(compound)

    ccl_df = pandas.read_table(metacelllinePath)
    ccl_df = ccl_df.set_index("master_ccl_id")

    response_df = pandas.read_table(responsePath)

    metaexperiment_df = pandas.read_table(metaexperimentPath)
    metaexperiment_df = metaexperiment_df.drop_duplicates("experiment_id")
    metaexperiment_df = metaexperiment_df.set_index("experiment_id")

    curve_df = pandas.read_table(curvePath)
    curve_df = curve_df.set_index(["experiment_id", "master_cpd_id"])

    project_compounds = {}
    for i, row in response_df.iterrows():
        exp_id = row['experiment_id']
        cpd_id = row['master_cpd_id']
        ccl_id = metaexperiment_df.loc[exp_id]['master_ccl_id']

        # remap sample name if possible
        ccl_name = ccl_df.loc[ccl_id]['ccl_name']
        if ccl_name in celllines:
            ccl_name = celllines[ccl_name]

        # Track drugs for project
        project_id = "CTRP_Unkown"
        if ccl_name in projects:
            project_id = "CTRP_{}".format(projects[ccl_name])
        proj = Project(project_id)
        if proj.gid() not in project_compounds:
            project_compounds[proj.gid()] = {}

        cpd_name = compound_df.loc[cpd_id]['cpd_name']
        auc = row['area_under_curve']
        ec50 = row['apparent_ec50_umol']

        curve_sub = curve_df.loc[(exp_id, cpd_id)]
        conc = curve_sub["cpd_conc_umol"]
        resp = curve_sub["cpd_avg_pv"]

        dr = DrugResponse(submitter_id=ccl_name, submitter_compound_id=cpd_name, source="ctrp",
                          act_area=auc, ec50=ec50, doses_um=list(conc),
                          activity_data_median=list(resp))
        emitter.emit_vertex(dr)
        compound = compound_factory(name=cpd_name)
        emitter.emit_edge(
            ResponseIn(),
            dr.gid(),
            Aliquot.make_gid("%s:DrugResponse:%s" % (ccl_name, cpd_name))
        )
        emitter.emit_edge(
            ResponseTo(),
            dr.gid(),
            compound.gid()
        )
        if compound.gid() not in project_compounds[proj.gid()]:
            emitter.emit_edge(
                TestedIn(),
                compound.gid(),
                proj.gid())
            project_compounds[proj.gid()][compound.gid()] = True

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()