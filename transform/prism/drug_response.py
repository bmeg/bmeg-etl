import json
import pandas
import sys

from bmeg.emitter import JSONEmitter
from bmeg import (Aliquot, DrugResponse, Project, Compound,
                  Compound_Projects_Project,
                  DrugResponse_Aliquot_Aliquot,
                  DrugResponse_Compounds_Compound)


def process_screen(screen_name, ex_df, fc_df, dr_df, compound_map, emitter):
    drugs = list(set(ex_df.broad_id))
    for drug_id in drugs:
        ex = ex_df[ex_df.broad_id == drug_id]
        doses = ex.dose.tolist()
        drug_name = ex.name[0]
        compound = compound_map.get(drug_name)
        if not compound:
            print("WARNING: compound not in lookup: {}".format(drug_name), file=sys.stderr)
            continue

        for cellline_id, row in fc_df.iterrows():
            if "FAILED" in cellline_id:
                continue

            # TODO: convert nan to None
            responses = row[ex.index].tolist()

            if isinstance(dr_df, pandas.DataFrame):
                try:
                    dr_vals = dr_df[(dr_df.broad_id == drug_id) & (dr_df.depmap_id == cellline_id)].iloc[0, ].dropna().to_dict()
                except Exception:
                    dr_vals = {}
            else:
                dr_vals = {}

            # create drug response vertex
            dr = DrugResponse(
                id=DrugResponse.make_gid("PRISM", screen_name, cellline_id, drug_name),
                einf=dr_vals.get("upper_limit"),
                ec50=dr_vals.get("ec50"),
                ic50=dr_vals.get("ic50"),
                aac=dr_vals.get("aac"),
                hs=dr_vals.get("slope"),
                dose_um=doses,
                response=responses,
                source_cell_name=cellline_id,
                source_drug_name=drug_name,
                project_id=Project.make_gid("PRISM")
            )
            emitter.emit_vertex(dr)

            emitter.emit_edge(
                DrugResponse_Aliquot_Aliquot(
                    from_gid=dr.gid(),
                    to_gid=Aliquot.make_gid("PRISM:%s" % (cellline_id))
                ),
                emit_backref=True
            )

            emitter.emit_edge(
                DrugResponse_Compounds_Compound(
                    from_gid=dr.gid(),
                    to_gid=compound.gid()
                ),
                emit_backref=True
            )

"""
def transform_compounds(drug_lookup_path='source/prism/compound_lookup.tsv',
                        emitter_prefix=None,
                        emitter_directory='prism'):

    emitter = JSONEmitter(prefix=emitter_prefix, directory=emitter_directory)
    with open(drug_lookup_path) as handle:
        for line in handle:
            row = line.split("\t")
            if row[1] == "":
                compound = Compound(
                    id=Compound.make_gid(drug_name),
                    project_id=Project.make_gid("Reference"),
                    id_source="prism", **json.loads(row[1])
                )
                compound['id'] = Compound.make_gid(compound['id'])
            else:
                compound = Compound(id=Compound.make_gid("NO_ONTOLOGY:%s" % (row[0])),
                                    id_source="NO_ONTOLOGY",
                                    submitter_id=row[0],
                                    project_id=Project.make_gid('Reference'))

            emitter.emit_vertex(compound)
            emitter.emit_edge(
                Compound_Projects_Project(
                    from_gid=compound.gid(),
                    to_gid=Project.make_gid("PRISM")
                ),
                emit_backref=True
            )
    emitter.close()
"""


def transform_compounds(
                        primary_file="source/prism/primary-screen-replicate-collapsed-treatment-info.csv",
                        secondary_file="source/prism/secondary-screen-replicate-collapsed-treatment-info.csv",
                        emitter_prefix=None,
                        emitter_directory='prism'):
    emitter = JSONEmitter(prefix=emitter_prefix, directory=emitter_directory)
    df1 = pandas.read_csv(primary_file)
    df2 = pandas.read_csv(secondary_file)
    for drug in set(df1.name.tolist() + df2.name.tolist()):
        if not pandas.isna(drug):
            compound = Compound(id=Compound.make_gid(drug),
                                id_source="prism",
                                submitter_id=drug,
                                project_id=Project.make_gid('Reference'))
            emitter.emit_vertex(compound)
            emitter.emit_edge(
                Compound_Projects_Project(
                    from_gid=compound.gid(),
                    to_gid=Project.make_gid("PRISM")
                ),
                emit_backref=True
            )
    emitter.close()



def transform_primary(drug_lookup_path='source/prism/compound_lookup.tsv',
                      primary_foldchange_path='source/prism/primary-screen-replicate-collapsed-logfold-change.csv',
                      primary_treatment_path='source/prism/primary-screen-replicate-collapsed-treatment-info.csv',
                      emitter_prefix="primary_screen",
                      emitter_directory='prism'):

    emitter = JSONEmitter(prefix=emitter_prefix, directory=emitter_directory)

    compound_map = {}
    with open(drug_lookup_path) as handle:
        for line in handle:
            row = line.split("\t")
            if row[1] == "":
                compound = Compound(**json.loads(row[1]))
                compound['id'] = Compound.make_gid(compound['id'])
            else:
                compound = Compound(id=Compound.make_gid(row[0]),
                                    id_source="prism",
                                    submitter_id=row[0],
                                    project_id=Project.make_gid('Reference'))

            compound_map[row[0]] = compound

    primary_fc = pandas.read_csv(primary_foldchange_path, index_col=0)
    primary_ex = pandas.read_csv(primary_treatment_path, index_col=0)

    process_screen("PRIMARY", primary_ex, primary_fc, None, compound_map, emitter)

    emitter.close()


def transform_secondary(drug_lookup_path='source/prism/compound_lookup.tsv',
                        secondary_foldchange_path='source/prism/secondary-screen-replicate-collapsed-logfold-change.csv',
                        secondary_response_path='source/prism/secondary-screen-dose-response-curve-parameters.csv',
                        secondary_treatment_path='source/prism/secondary-screen-replicate-collapsed-treatment-info.csv',
                        emitter_prefix="secondary_screen",
                        emitter_directory='prism'):

    emitter = JSONEmitter(prefix=emitter_prefix, directory=emitter_directory)

    compound_map = {}
    with open(drug_lookup_path) as handle:
        for line in handle:
            row = line.split("\t")
            if row[1] == "":
                compound = Compound(**json.loads(row[1]))
                compound['id'] = Compound.make_gid(compound['id'])
            else:
                compound = Compound(id=Compound.make_gid(row[0]),
                                    id_source="prism",
                                    submitter_id=row[0],
                                    project_id=Project.make_gid('Reference'))

            compound_map[row[0]] = compound

    secondary_dr = pandas.read_csv(secondary_response_path, low_memory=False).astype(object)
    secondary_dr["aac"] = 1 + -secondary_dr.auc
    secondary_fc = pandas.read_csv(secondary_foldchange_path, index_col=0).astype(object)
    secondary_ex = pandas.read_csv(secondary_treatment_path, index_col=0).astype(object)

    process_screen("SECONDARY", secondary_ex, secondary_fc, secondary_dr, compound_map, emitter)

    emitter.close()
