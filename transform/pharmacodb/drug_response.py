import pandas
import bmeg.ioutils

from bmeg.emitter import JSONEmitter
from bmeg import (Aliquot, DrugResponse, Project,
                  Compound_Projects_Project,
                  DrugResponse_Aliquot_Aliquot,
                  DrugResponse_Compounds_Compound)
from bmeg.enrichers.drug_enricher import compound_factory


def transform(cellline_lookup_path='source/ccle/cellline_id_lookup.tsv',
              cells_path='source/pharmacodb/cells.tsv.gz',
              drugs_path='source/pharmacodb/drugs.tsv.gz',
              drug_annots_path='source/pharmacodb/drug_annots.tsv.gz',
              experiments_path='source/pharmacodb/experiments.tsv.gz',
              dose_response_path='source/pharmacodb/dose_responses.tsv.gz',
              profiles_path='source/pharmacodb/profiles.tsv.gz',
              emitter_prefix=None,
              emitter_directory='pharmacodb'):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    cellline_id_lookup = bmeg.ioutils.read_lookup(cellline_lookup_path)

    dose_response = pandas.read_csv(dose_response_path, sep="\t")

    cells = pandas.read_csv(cells_path, sep="\t")
    drugs = pandas.read_csv(drugs_path, sep="\t")
    drug_annots = pandas.read_csv(drug_annots_path, sep="\t")
    experiments = pandas.read_csv(experiments_path, sep="\t")
    profiles = pandas.read_csv(profiles_path, sep="\t")

    merged_data = pandas.merge(
        cells[['cell_id', 'cell_name']],
        pandas.merge(
            profiles,
            pandas.merge(
                experiments,
                pandas.merge(
                    drugs[['drug_id', 'drug_name']],
                    drug_annots,
                    on='drug_id'
                ),
                on='drug_id'
            ),
            on='experiment_id'
        ),
        on='cell_id'
    )

    # dataset_id,dataset_name
    # 1,CCLE
    # 2,CTRPv2
    # 5,GDSC1000
    merged_data = merged_data[merged_data.dataset_id.isin([1, 2, 5])]
    merged_data['dataset_name'] = merged_data.dataset_id.map({1: 'CCLE', 2: 'CTRP', 5: 'GDSC'})

    emitted_compounds = {}
    for index, row in merged_data.iterrows():
        row = row.dropna().to_dict()
        cell_name = cellline_id_lookup.get(row["cell_name"], row["cell_name"])
        dr = DrugResponse(
            id=DrugResponse.make_gid(row["dataset_name"], row["experiment_id"],
                                     cell_name, row["drug_name"]),
            einf=row.get("Einf"),
            ec50=row.get("EC50"),
            ic50=row.get("IC50"),
            aac=row.get("AAC"),
            hs=row.get("HS"),
            dss1=row.get("DSS1"),
            dss2=row.get("DSS2"),
            dss3=row.get("DSS3"),
            dose_um=dose_response[dose_response.experiment_id == row["experiment_id"]].dose.tolist(),
            response=dose_response[dose_response.experiment_id == row["experiment_id"]].response.tolist(),
            source_cell_name=row.get("cell_name"),
            source_drug_name=row.get("drug_name"),
            project_id=Project.make_gid(row["dataset_name"])
        )
        emitter.emit_vertex(dr)

        drug_name = row["drug_name"]
        if row.get("pubchem"):
            drug_name = "CID{}".format(int(row.get("pubchem")))
        compound = compound_factory(name=drug_name)

        if compound.gid() not in emitted_compounds:
            emitter.emit_vertex(compound)
            emitter.emit_edge(
                Compound_Projects_Project(
                    from_gid=compound.gid(),
                    to_gid=Project.make_gid(row["dataset_name"])
                ),
                emit_backref=True
            )
            emitted_compounds[compound.gid()] = None

        # add edge to compound
        emitter.emit_edge(
            DrugResponse_Compounds_Compound(
                from_gid=dr.gid(),
                to_gid=compound.gid()
            ),
            emit_backref=True,
        )
        # add edge to aliquot
        emitter.emit_edge(
            DrugResponse_Aliquot_Aliquot(
                from_gid=dr.gid(),
                to_gid=Aliquot.make_gid("{}:{}".format(row["dataset_name"], cell_name))
            ),
            emit_backref=True
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
