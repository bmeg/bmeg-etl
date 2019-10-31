import pandas
import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import (Aliquot, DrugResponse, Project,
                  Compound_Projects_Project,
                  DrugResponse_Aliquot_Aliquot,
                  DrugResponse_Compounds_Compound)
from bmeg.enrichers.drug_enricher import compound_factory


def transform(cellline_lookup_path='source/ccle/cellline_lookup.tsv',
              cellosaurus_path='source/pharmacodb/cellosaurus.csv',
              cells_path='source/pharmacodb/cells.csv',
              drugs_path='source/pharmacodb/drugs.csv',
              drug_annots_path='source/pharmacodb/drug_annots.csv',
              experiments_path='source/pharmacodb/experiments.csv',
              dose_response_path='source/pharmacodb/dose_responses.csv',
              profiles_path='source/pharmacodb/profiles.csv',
              source_cell_path='source/pharmacodb/source_cell_names.csv',
              source_drug_path='source/pharmacodb/source_drug_names.csv',
              emitter_prefix=None,
              emitter_directory='pharmacodb'):

    cellline_lookup_path='source/ccle/cellline_lookup.tsv'
    cellosaurus_path='source/pharmacodb/cellosaurus.csv'
    cells_path='source/pharmacodb/cells.csv'
    drugs_path='source/pharmacodb/drugs.csv'
    drug_annots_path='source/pharmacodb/drug_annots.csv'
    experiments_path='source/pharmacodb/experiments.csv'
    dose_response_path='source/pharmacodb/dose_responses.csv'
    profiles_path='source/pharmacodb/profiles.csv'
    source_cell_path='source/pharmacodb/source_cell_names.csv'
    source_drug_path='source/pharmacodb/source_drug_names.csv'

    celllines = bmeg.ioutils.read_lookup(cellline_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    dose_response = pandas.read_csv(dose_response_path)

    cellosaurus = pandas.read_csv(cellosaurus_path)
    cells = pandas.read_csv(cells_path)
    drugs = pandas.read_csv(drugs_path)
    drug_annots = pandas.read_csv(drug_annots_path)
    experiments = pandas.read_csv(experiments_path)
    profiles = pandas.read_csv(profiles_path)
    source_cell_names = pandas.read_csv(source_cell_path)
    source_drug_names = pandas.read_csv(source_drug_path)

    cell_names = pandas.concat([
        pandas.merge(
            cells[["cell_id", "cell_name"]],
            cellosaurus.rename(columns={"identifier": "cell_name"})[["cell_name", "sy", "dr"]],
            on='cell_name',
            how='inner'
        ),
        pandas.merge(
            source_cell_names[["cell_id", "cell_name"]],
            cellosaurus.rename(columns={"identifier": "cell_name"})[["cell_name", "sy", "dr"]],
            on='cell_name',
            how='inner'
        )
    ]).drop_duplicates()

    def getSynonyms(cell_table, cell_id):
        entry = cell_table[cell_table.cell_id == cell_id]
        if entry.shape[0] == 0:
            return []
        syn = []
        syn.append(entry.cell_name.tolist()[0])
        try:
            syn = syn + entry.sy.map(lambda s: [x.strip() for x in str(s).split(";")]).tolist()[0]
        except Exception:
            pass
        try:
            syn = syn + entry.dr.map(lambda s: [x.split(";")[1].strip() for x in str(s).split("|")]).tolist()[0]
        except Exception:
            pass
        return list(set(syn))

    def lookupBroadId(name_aliases, broad_lookup, name, cell_id):
        names = [name, str(name), str(name).lower(), str(name).upper()]
        name_check = [n in broad_lookup for n in names]
        if any(name_check):
            return broad_lookup[names[name_check.index(True)]]

        else:
            syns = getSynonyms(name_aliases, cell_id)
            for s in syns:
                if s in celllines:
                    return broad_lookup[s]

            name_stripped = str(name).replace("-", "").replace(".", "").replace("/", "").replace(" ", "")
            stripped_names = [name_stripped, name_stripped.lower(), name_stripped.upper()]
            stripped_name_check = [n in broad_lookup for n in stripped_names]
            if any(stripped_name_check):
                return broad_lookup[stripped_names[stripped_name_check.index(True)]]

            else:
                return None

    merged_data = pandas.merge(
        cells[["cell_id", "cell_name"]],
        pandas.merge(
            profiles,
            pandas.merge(
                experiments,
                pandas.merge(
                    source_drug_names.rename(columns={"drug_name": "source_drug_name"}, inplace=False)[["drug_id", "source_drug_name"]],
                    pandas.merge(
                        drugs[["drug_id", "drug_name"]],
                        drug_annots,
                        on="drug_id"
                    ),
                    on="drug_id"
                ),
                on="drug_id"
            ),
            on="experiment_id"
        ),
        on="cell_id"
    )

    # dataset_id,dataset_name
    # 1,CCLE
    # 2,CTRPv2
    # 5,GDSC1000
    merged_data = merged_data[merged_data.dataset_id.isin([1, 2, 5])]
    merged_data["dataset_name"] = merged_data.dataset_id.map({1: "CCLE", 2: "CTRP", 5: "GDSC"})

    found = 0
    missing = 0
    missing_names = []
    for index, row in merged_data[["cell_id", "cell_name"]].drop_duplicates().iterrows():
        name = lookupBroadId(cell_names, celllines, row["cell_name"], row["cell_id"])
        if name:
            found += 1
        else:
            missing += 1
            missing_names.append(row["cell_name"])

    # print("found:", found)
    # 1275
    # print("missing:", missing)
    # 104
    # merged_data[merged_data.cell_name.isin(missing_names)][["cell_name", "dataset_name"]].drop_duplicates()["dataset_name"].value_counts()
    # GDSC    57
    # CTRP    46
    # CCLE     4

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
