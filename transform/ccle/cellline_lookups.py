import os
from pydash import is_blank

import bmeg.ioutils


def create_project_lookup(path="source/ccle/DepMap-2019q1-celllines.csv_v2.csv",
                          outdir="source/ccle/"):
    """
    Create lookup for DepMap cell line to project code

    Lookup key is DepMap_ID
    """
    lookup = {}
    input_stream = bmeg.ioutils.read_csv(path)
    for line in input_stream:
        if "NORMAL" in line["CCLE_Name"] or "MERGED" in line["CCLE_Name"]:
            continue

        # TODO: convert to TCGA project short codes
        project_id = "_".join(line["Primary Disease"].split())
        # project_id = line["Subtype Disease"]

        lookup[line["DepMap_ID"]] = project_id

    assert len(lookup.keys()), "Project lookup is empty!"

    handle = open(os.path.join(outdir, "cellline_project_lookup.tsv"), 'w')
    for key, value in lookup.items():
        handle.write("{}\t{}\n".format(key, value))
    handle.close()
    return


def create_phenotype_lookup(path="source/ccle/DepMap-2019q1-celllines.csv_v2.csv",
                            outdir="source/ccle/"):
    """
    Create lookup for DepMap cell line to phenotype

    Lookup key is DepMap_ID
    """
    lookup = {}
    input_stream = bmeg.ioutils.read_csv(path)
    for line in input_stream:
        if "NORMAL" in line["CCLE_Name"] or "MERGED" in line["CCLE_Name"]:
            continue

        phenotype_name = line.get('Subtype Disease', None)
        if not phenotype_name or is_blank(phenotype_name):
            phenotype_name = line.get('Primary Disease')
        lookup[line["DepMap_ID"]] = phenotype_name

    assert len(lookup.keys()), "PHenotype lookup is empty!"

    handle = open(os.path.join(outdir, "cellline_phenotype_lookup.tsv"), 'w')
    for key, value in lookup.items():
        handle.write("{}\t{}\n".format(key, value))
    handle.close()
    return


def create_cellline_lookup(path="source/ccle/DepMap-2019q1-celllines.csv_v2.csv",
                           outdir="source/ccle/"):
    """
    create lookup of cell line names to DepMap_ID
    """

    lookup = {}
    input_stream = bmeg.ioutils.read_csv(path)
    for line in input_stream:
        if "NORMAL" in line["CCLE_Name"] or "MERGED" in line["CCLE_Name"]:
            continue

        lookup[line["DepMap_ID"]] = line["DepMap_ID"]
        lookup[line["CCLE_Name"]] = line["DepMap_ID"]
        lookup[line["COSMIC_ID"]] = line["DepMap_ID"]
        lookup[line["Sanger ID"]] = line["DepMap_ID"]

        aliases = []

        prefix = line["CCLE_Name"].split('_')[0]
        blacklisted_prefixes = ["KMH2", "MS1", "TT"]
        if prefix not in blacklisted_prefixes:
            aliases.append(prefix)

        split_aliases = line["Aliases"].split(";")
        if len(split_aliases) > 1:
            for x in split_aliases:
                aliases.append(x)

        for a in aliases:
            lookup[a] = line["DepMap_ID"]

    assert len(lookup.keys()), "Cell line lookup is empty!"

    handle = open(os.path.join(outdir, "cellline_lookup.tsv"), 'w')
    for key, value in lookup.items():
        handle.write("{}\t{}\n".format(key, value))
    handle.close()
    return


if __name__ == "__main__":
    create_cellline_lookup()
    create_phenotype_lookup()
    create_project_lookup()
