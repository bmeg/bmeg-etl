import pandas
from bmeg.vertex import DrugResponse, Aliquot, Case, Project, Program
from bmeg.edge import ResponseIn, ResponseTo, InProject, InProgram, TestedIn
from bmeg.enrichers.drug_enricher import compound_factory
from bmeg.emitter import JSONEmitter
from bmeg.ccle import build_ccle2depmap_conversion_table, build_project_lookup, missing_ccle_cellline_factory


def transform(sample_path='outputs/ccle/Sample.Vertex.json.gz',
              metadrugPath='source/ctrp/v20.meta.per_compound.txt',
              metacelllinePath="source/ctrp/v20.meta.per_cell_line.txt",
              responsePath="source/ctrp/v20.data.curves_post_qc.txt",
              metaexperimentPath="source/ctrp/v20.meta.per_experiment.txt",
              curvePath="source/ctrp/v20.data.per_cpd_post_qc.txt",
              emitter_prefix='ctrp',
              emitter_directory='ctrp'):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    prog = Program(program_id="CTRP")
    emitter.emit_vertex(prog)

    # lookup table to convert CCLE names to DepMap_IDs
    ccle_id_lookup = build_ccle2depmap_conversion_table(sample_path)
    # lookup table for projects
    projects = build_project_lookup(sample_path)

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

    missing_cell_lines = []
    case_gids = []
    project_gids = []
    project_compounds = {}
    for i, row in response_df.iterrows():
        exp_id = row['experiment_id']
        cpd_id = row['master_cpd_id']
        ccl_id = metaexperiment_df.loc[exp_id]['master_ccl_id']

        # remap sample name if possible
        ccl_name = ccl_df.loc[ccl_id]['ccl_name']
        if ccl_name in ccle_id_lookup:
            ccl_name = ccle_id_lookup[ccl_name]
        else:
            if ccl_name not in missing_cell_lines:
                missing_cell_lines.append(ccl_name)

        # Create project and link to case and program
        project_id = "CTRP_Unkown"
        if ccl_name in projects:
            project_id = "CTRP_{}".format(projects[ccl_name])
        proj = Project(project_id)
        if proj.gid() not in project_gids:
            emitter.emit_vertex(proj)
            emitter.emit_edge(
                InProgram(),
                proj.gid(),
                prog.gid()
            )
            project_gids.append(proj.gid())
            project_compounds[proj.gid()] = {}
        if Case.make_gid(ccl_name) not in case_gids and ccl_name not in missing_cell_lines:
            emitter.emit_edge(
                InProject(),
                Case.make_gid(ccl_name),
                proj.gid(),
            )
            case_gids.append(Case.make_gid(ccl_name))

        cpd_name = compound_df.loc[cpd_id]['cpd_name']
        auc = row['area_under_curve']
        ec50 = row['apparent_ec50_umol']

        # curve_sub = curve_df.loc[ (curve_df["experiment_id"] == exp_id) & (curve_df["master_cpd_id"] == cpd_id) ]
        curve_sub = curve_df.loc[(exp_id, cpd_id)]
        conc = curve_sub["cpd_conc_umol"]
        resp = curve_sub["cpd_avg_pv"]

        dr = DrugResponse(sample_id=ccl_name, compound_id=cpd_name, source="CTRP",
                          act_area=auc, ec50=ec50, doses_um=list(conc),
                          activity_data_median=list(resp))
        emitter.emit_vertex(dr)
        compound = compound_factory(name=cpd_name)
        emitter.emit_edge(
            ResponseIn(),
            dr.gid(),
            Aliquot.make_gid(ccl_name)
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

    # generate project, case, sample, aliquot for missing cell lines
    missing_ccle_cellline_factory(emitter=emitter,
                                  missing_ids=missing_cell_lines,
                                  project_id="CTRP_Unknown",
                                  project_gids=project_gids)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
