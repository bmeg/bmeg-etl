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
              drugs_path='source/pharmacodb/drugs.csv',
              drug_annots_path='source/pharmacodb/drug_annots.csv',
              experiments_path='source/pharmacodb/experiments.csv',
              dose_response_path='source/pharmacodb/dose_responses.csv',
              profiles_path='source/pharmacodb/profiles.csv',
              source_drug_path='source/pharmacodb/source_drug_names.csv',
              emitter_prefix=None,
              emitter_directory='pharmacodb'):

    cellline_lookup_path='source/ccle/cellline_lookup.tsv'

    cells_path='source/pharmacodb/cells.csv'
    drugs_path='source/pharmacodb/drugs.csv'
    drug_annots_path='source/pharmacodb/drug_annots.csv'
    experiments_path='source/pharmacodb/experiments.csv'
    dose_response_path='source/pharmacodb/dose_responses.csv'
    profiles_path='source/pharmacodb/profiles.csv'
    source_drug_path='source/pharmacodb/source_drug_names.csv'

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    dose_response = pandas.read_csv(dose_response_path)

    cells = pandas.read_csv(cells_path)
    drugs = pandas.read_csv(drugs_path)
    drug_annots = pandas.read_csv(drug_annots_path)
    experiments = pandas.read_csv(experiments_path)
    profiles = pandas.read_csv(profiles_path)
    source_drug_names = pandas.read_csv(source_drug_path)

    merged_data = pandas.merge(
        cells[['cell_id', 'cell_name']],
        pandas.merge(
            profiles,
            pandas.merge(
                experiments,
                pandas.merge(
                    source_drug_names.rename(columns={'drug_name': 'source_drug_name'}, inplace=False)[['drug_id', 'source_drug_name']],
                    pandas.merge(
                        drugs[['drug_id', 'drug_name']],
                        drug_annots,
                        on='drug_id'
                    ),
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

    for index, row in merged_data.iterrows():
        cell_name = celllines.get(row.cell_name, row.cell_name)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
