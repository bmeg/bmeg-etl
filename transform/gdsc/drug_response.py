import pandas

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Aliquot, DrugResponse, Project
from bmeg.edge import ResponseIn, ResponseTo, TestedIn
from bmeg.enrichers.drug_enricher import compound_factory


def transform(project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              drugs_meta_path="source/gdsc/Screened_Compounds.xlsx",
              ic50_path="source/gdsc/GDSC_IC50.csv",
              auc_path="source/gdsc/GDSC_AUC.csv",
              emitter_prefix="gdsc",
              emitter_directory="gdsc"):

    projects = bmeg.ioutils.read_lookup(project_lookup_path)

    emitter = JSONEmitter(prefix=emitter_prefix, directory=emitter_directory)

    drugs_df = pandas.read_excel(drugs_meta_path, sheet_name=0, header=0)
    drugs_df.set_index("DRUG_ID", inplace=True, drop=True)

    auc_df = pandas.read_csv(auc_path, sep=",", index_col=0)
    auc_df.index = [int(x.replace('GDSC:', '')) for x in auc_df.index.tolist()]

    ic50_df = pandas.read_csv(ic50_path, sep=",", index_col=0)
    ic50_df.index = [int(x.replace('GDSC:', '')) for x in ic50_df.index.tolist()]

    emitted_compounds = {}
    project_compounds = {}
    for drug_id, row in ic50_df.iterrows():
        drug_name = drugs_df.loc[drug_id]['DRUG_NAME']
        for cellline_id, ic50_val in row.to_dict().items():
            auc_val = auc_df.loc[drug_id][cellline_id]

            # Track drugs for project
            project_id = "GDSC_%s" % (projects.get(cellline_id, "Unknown"))
            proj = Project(project_id)
            if proj.gid() not in project_compounds:
                project_compounds[proj.gid()] = {}

            dr = DrugResponse(submitter_id=cellline_id,
                              submitter_compound_id=drug_name,
                              source="GDSC",
                              act_area=auc_val,
                              ic50=ic50_val)
            emitter.emit_vertex(dr)
            compound = compound_factory(name=drug_name)
            if compound.gid() not in emitted_compounds:
                emitter.emit_vertex(compound)
                emitted_compounds[compound.gid()] = True
            emitter.emit_edge(
                ResponseIn(),
                dr.gid(),
                Aliquot.make_gid("%s:DrugResponse:%s" % (cellline_id, drug_name))
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


if __name__ == "__main__":
    transform()
