import json
import os
import pandas
import re
from pydash import is_blank


def create_cellline_lookups(depmap_path="source/ccle/DepMap-2019q1-celllines.csv_v2.csv",
                            cellosaurus_path="source/pharmacodb/cellosaurus.csv",
                            cells_path="source/pharmacodb/cells.csv",
                            source_cell_path="source/pharmacodb/source_cell_names.csv",
                            outdir="source/ccle/"):
    """
    create cell line lookup tables:
      - map cell name aliases to DepMap IDs
      - map DepMap IDs to phenotypes
    """

    lookup = {}
    depmap = pandas.read_csv(depmap_path, dtype=str)
    for index, line in depmap.iterrows():
        line = line.dropna()

        if "MERGED" in line["CCLE_Name"]:
            m = re.search("\[MERGED_TO_(ACH-[0-9]+)\](.*)", line["CCLE_Name"])
            lookup[line["DepMap_ID"]] = m.group(1)
            lookup[m.group(2)] = m.group(1)
            if line.get("COSMIC_ID"):
                lookup[line["COSMIC_ID"]] = m.group(1)
            continue

        keys = ["DepMap_ID", "CCLE_Name", "COSMIC_ID"]
        for k in keys:
            if line.get(k) and not is_blank(line.get(k)):
                lookup[line.get(k)] = line["DepMap_ID"]

        if "MATCHED_NORMAL" in line["CCLE_Name"]:
            continue

        aliases = []
        prefix = line["CCLE_Name"].split("_")[0]
        blacklisted_prefixes = ["KMH2", "MS1", "TT"]
        if prefix not in blacklisted_prefixes:
            aliases.append(prefix)

        if line.get("Aliases") and not is_blank(line.get(k)):
            split_aliases = line.get("Aliases").split(";")
            for x in split_aliases:
                aliases.append(x)

        for a in aliases:
            lookup[a] = line["DepMap_ID"]

    assert len(lookup.keys()), "Cell line lookup is empty!"

    # map pharmacodb cell names to Broad IDs
    cellosaurus = pandas.read_csv(cellosaurus_path)
    cells = pandas.read_csv(cells_path)
    source_cell_names = pandas.read_csv(source_cell_path)
    pharmacodb_aliases = pandas.merge(
        pandas.merge(
            cells[["cell_id", "cell_name"]],
            source_cell_names.groupby("cell_id")["cell_name"].apply(lambda x: list(set(x))).reset_index(name="project_aliases"),
            on="cell_id",
            how="outer"
        ),
        cellosaurus.rename(columns={"identifier": "cell_name"})[["cell_name", "sy", "dr", "di"]],
        on="cell_name",
        how="left"
    )
    pharmacodb_aliases["project_aliases"] = pharmacodb_aliases.project_aliases.apply(lambda x: x if x == x else [])
    pharmacodb_aliases["sy"] = pharmacodb_aliases.sy.apply(lambda s: [x.strip() for x in str(s).split(";") if x != "nan"])
    pharmacodb_aliases["dr"] = pharmacodb_aliases.dr.apply(lambda s: [x.split(';')[1].strip() for x in str(s).split('|') if ';' in x])
    pharmacodb_aliases["ncit_id"] = pharmacodb_aliases.di.apply(lambda x: x.split(";")[1] if isinstance(x, str) else x)
    pharmacodb_aliases["ncit_desc"] = pharmacodb_aliases.di.apply(lambda x: x.split(";")[2] if isinstance(x, str) else x)

    for index, row in pharmacodb_aliases.iterrows():
        name = row.cell_name
        names = []
        names = [name, str(name), str(name).lower(), str(name).upper()]
        name_check = [n in lookup for n in names]
        if any(name_check):
            lookup[name] = lookup[names[name_check.index(True)]]
            continue

        try:
            syns = row.sy + row.dr + row.project_aliases
        except Exception as e:
            print(row)
            raise e
        for s in syns:
            if s in lookup:
                lookup[name] = lookup[s]
                continue

        name_stripped = str(name).replace('-', '').replace('.', '').replace('/', '').replace(' ', '')
        stripped_names = [name_stripped, name_stripped.lower(), name_stripped.upper()]
        stripped_name_check = [n in lookup for n in stripped_names]
        if any(stripped_name_check):
            lookup[name] = lookup[stripped_names[stripped_name_check.index(True)]]
            continue

    handle = open(os.path.join(outdir, "cellline_id_lookup.tsv"), "w")
    for key, value in lookup.items():
        handle.write("{}\t{}\n".format(key, value))
    handle.close()

    # create phenotype lookup
    phenotypes = {}
    for index, row in pharmacodb_aliases.iterrows():
        n = row.cell_name
        if n in lookup:
            n = lookup[n]
        if not pandas.isna(row.ncit_desc):
            phenotypes[n] = row.ncit_desc

    for index, line in depmap.iterrows():
        broad_id = line["DepMap_ID"]
        if broad_id not in phenotypes:
            phenotype_name = line["Primary Disease"]
            phenotypes[broad_id] = phenotype_name.lower()

    assert len(phenotypes.keys()), "Phenotype lookup is empty!"

    handle = open(os.path.join(outdir, "cellline_phenotype_lookup.tsv"), "w")
    for key, value in phenotypes.items():
        handle.write("{}\t{}\n".format(key, value))
    handle.close()

    # create cell line properties lookup
    cell_props = {}
    for index, line in depmap.iterrows():
        broad_id = line["DepMap_ID"]
        props = line.to_dict()
        cell_props[broad_id] = json.dumps(props)

    assert len(cell_props.keys()), "Cell properties lookup is empty!"

    handle = open(os.path.join(outdir, "cellline_properties_lookup.tsv"), "w")
    for key, value in cell_props.items():
        handle.write("{}\t{}\n".format(key, value))
    handle.close()

    return


if __name__ == "__main__":
    create_cellline_lookups()
