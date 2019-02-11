""" match drug response data to aliquot """

from types import SimpleNamespace as SN

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import DrugResponse, Aliquot, Individual, Project, Program
from bmeg.edge import ResponseIn, ResponseTo, InProject, InProgram
from bmeg.enrichers.drug_enricher import compound_factory
from bmeg.ccle import build_ccle2depmap_conversion_table, build_project_lookup, missing_ccle_cellline_factory


NAMES = {
    'CCLE Cell Line Name': 'sample_id',
    # 'Primary Cell Line Name': 'primary_cell_line_name',
    'Compound': 'compound_id',
    # 'Target': 'target',
    'Doses (uM)': 'doses_um',
    'Activity Data (median)': 'activity_data_median',
    'Activity SD': 'activity_sd',
    'Num Data': 'num_data',
    'FitType': 'fit_type',
    'EC50 (uM)': 'ec50',
    'IC50 (uM)': 'ic50',
    'Amax': 'amax',
    'ActArea': 'act_area',
}


def transform(biosample_path='outputs/ccle/Biosample.Vertex.json.gz',
              drug_response_path='source/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv',
              emitter_prefix='drug_response',
              emitter_directory="ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    prog = Program(program_id="CCLE")
    emitter.emit_vertex(prog)

    # lookup table to convert CCLE names to DepMap_IDs
    samples = build_ccle2depmap_conversion_table(biosample_path)
    # lookup table for projects
    projects = build_project_lookup(biosample_path)

    # input and map
    input_stream = bmeg.ioutils.read_csv(drug_response_path)
    floats = ['amax', 'act_area', 'ec50', 'ic50', 'num_data']
    compound_gids = []
    # read the drug response csv
    missing_cell_lines = []
    individual_gids = []
    project_gids = []
    for line in input_stream:
        # map the names to snake case
        mline = {"source": "CCLE"}
        for k, v in NAMES.items():
            # line[v] = line.pop(k)
            # csv field to float array
            if ',' in line[k]:
                mline[v] = [float(s) for s in line[k].split(',')]
            # NA -> None
            elif line[k] == 'NA':
                mline[v] = None
            # string scalar to float
            elif line[k] and v in floats:
                mline[v] = float(line[k])
            else:
                mline[v] = line[k]
        drug_response = SN(**mline)

        sample_id = drug_response.sample_id
        # if no match, we will need to create project->individual->biosample->aliquot
        if sample_id in samples:
            drug_response.sample_id = samples[sample_id]
        elif sample_id.split("_")[0] in samples:
            drug_response.sample_id = samples[sample_id.split("_")[0]]
        else:
            if sample_id not in missing_cell_lines:
                missing_cell_lines.append(sample_id)

        # Create project and link to individual and program
        project_id = "CCLE_Unkown"
        if sample_id in projects:
            project_id = "CCLE_{}".format(projects[sample_id])
        elif "_".join(sample_id.split("_")[1:]) in projects:
            project_id = "CCLE_{}".format(projects["_".join(sample_id.split("_")[1:])])
        proj = Project(project_id)
        if Individual.make_gid(sample_id) not in individual_gids:
            emitter.emit_edge(
                InProject(),
                Individual.make_gid(sample_id),
                proj.gid(),
            )
            individual_gids.append(Individual.make_gid(sample_id))
        if proj.gid() not in project_gids:
            emitter.emit_vertex(proj)
            emitter.emit_edge(
                InProgram(),
                proj.gid(),
                prog.gid()
            )
            project_gids.append(proj.gid())

        # create drug_response vertex
        drug_resp = DrugResponse(**drug_response.__dict__)
        emitter.emit_vertex(drug_resp)
        #  and edge to aliquot
        emitter.emit_edge(
            ResponseIn(),
            drug_resp.gid(),
            Aliquot.make_gid(drug_response.sample_id),
        )
        # create compound
        compound = compound_factory(name=drug_response.compound_id)
        if compound.gid() not in compound_gids:
            emitter.emit_vertex(compound)
            compound_gids.append(compound.gid())
        # and an edge to it
        emitter.emit_edge(
            ResponseTo(),
            drug_resp.gid(),
            compound.gid(),
        )

    # generate project, individual, biosample, aliquot for missing cell lines
    missing_ccle_cellline_factory(emitter=emitter,
                                  missing_ids=missing_cell_lines,
                                  project_prefix="CCLE",
                                  individual_gids=individual_gids,
                                  project_gids=project_gids)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
