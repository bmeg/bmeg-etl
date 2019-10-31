import pandas
import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import (Aliquot, DrugResponse, Project,
                  Compound_Projects_Project,
                  DrugResponse_Aliquot_Aliquot,
                  DrugResponse_Compounds_Compound)
from bmeg.enrichers.drug_enricher import compound_factory


def transform(cellline_lookup_path='source/ccle/cellline_lookup.tsv',
              cells_path='source/pharmacodb/cells.csv',
              drug_annots_path='source/pharmacodb/drug_annots.csv',
              experiments_path='source/pharmacodb/experiments.csv',
              dose_response_path='source/pharmacodb/dose_responses.csv',
              profiles_path='source/pharmacodb/profiles.csv',
              source_cell_path='source/pharmacodb/source_cell_names.csv',
              source_drug_path='source/pharmacodb/source_drug_names.csv',
              emitter_prefix=None,
              emitter_directory='pharmacodb'):

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)
    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    cells = pandas.read_csv(cells_path)
    drug_annots = pandas.read_csv(drug_annots_path)
    experiments = pandas.read_csv(experiments_path)
    dose_response = pandas.read_csv(dose_response_path)
    profiles = pandas.read_csv(profiles_path)
    source_cell_names = pandas.read_csv(source_cell_path)
    source_drug_names = pandas.read_csv(source_drug_path)

    merged_data = pandas.merge(
        source_cell_names.rename(columns={"cell_name": "source_cell_name"}, inplace=False)[["cell_id", "source_cell_name"]],
        pandas.merge(
            cells[["cell_id", "cell_name"]],
            pandas.merge(
                profiles,
                pandas.merge(
                    experiments,
                    pandas.merge(
                        source_drug_names[["drug_id", "drug_name"]],
                        drug_annots,
                        on="drug_id"
                    ),
                    on="drug_id"
                ),
                on="experiment_id"
            ),
            on="cell_id"
        ),
        on="cell_id"
    )

    # dataset_id,dataset_name
    # 1,CCLE
    # 2,CTRPv2
    # 5,GDSC1000
    merged_data = merged_data[merged_data.dataset_id.isin([1, 2, 5])]

    names = merged_data.cell_name.unique().tolist()
    # found 633
    [name in celllines for name in names].count(True)
    # 746 missing
    [name in celllines for name in names].count(False)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
