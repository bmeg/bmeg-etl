import pandas

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import (Aliquot, DrugResponse, Project,
                  Compound_Projects_Project,
                  DrugResponse_Aliquot_Aliquot,
                  DrugResponse_Compounds_Compound)
from bmeg.enrichers.drug_enricher import compound_factory


def transform(cellline_lookup_path="source/ccle/cellline_lookup.tsv",
              project_lookup_path="source/ccle/cellline_project_lookup.tsv",
              drugs_meta_path="source/gdsc/Screened_Compounds.xlsx",
              ic50_path="source/gdsc/GDSC_IC50.csv",
              auc_path="source/gdsc/GDSC_AUC.csv",
              emitter_prefix="drug_response",
              emitter_directory="gdsc"):

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
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

            # correct for merged broad ids
            cellline_id = celllines.get(cellline_id, cellline_id)

            # Track drugs for project
            project_id = "GDSC_%s" % (projects.get(cellline_id, "Unknown"))
            proj = Project(submitter_id=Project.make_gid(project_id),
                           project_id=project_id)
            if proj.gid() not in project_compounds:
                project_compounds[proj.gid()] = {}

            # create drug response vertex
            dr = DrugResponse(submitter_id=DrugResponse.make_gid("GDSC", cellline_id, drug_name),
                              submitter_compound_id=drug_name,
                              auc=auc_val if not pandas.isnull(auc_val) else None,
                              ic50=ic50_val if not pandas.isnull(ic50_val) else None,
                              project_id=proj.gid())
            emitter.emit_vertex(dr)
            #  and edge to aliquot
            emitter.emit_edge(
                DrugResponse_Aliquot_Aliquot(
                    from_gid=dr.gid(),
                    to_gid=Aliquot.make_gid("GDSC:%s:DrugResponse:%s" % (cellline_id, drug_name))
                ),
                emit_backref=True
            )

            # create compound
            compound = compound_factory(name=drug_name)
            if compound.gid() not in emitted_compounds:
                emitter.emit_vertex(compound)
                emitted_compounds[compound.gid()] = True
            # and an edge to it
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


if __name__ == "__main__":
    transform()
