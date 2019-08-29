import os
import re
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
        # TODO: convert to TCGA project short codes
        project_id = "_".join(line["Primary Disease"].split())
        # project_id = line["Subtype Disease"]
        if not project_id or is_blank(project_id):
            continue
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
        phenotype_name = line.get('Primary Disease', None)
        if phenotype_name is None or is_blank(phenotype_name) or phenotype_name == "unknown":
            phenotype_name = "Unknown"

        lookup[line["DepMap_ID"]] = phenotype_name

    assert len(lookup.keys()), "Phenotype lookup is empty!"

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
        if "MERGED" in line["CCLE_Name"]:
            m = re.search("\[MERGED_TO_(ACH-[0-9]+)\](.*)", line["CCLE_Name"])
            lookup[line["DepMap_ID"]] = m.group(1)
            lookup[m.group(2)] = m.group(1)
            lookup[line["COSMIC_ID"]] = m.group(1)
            continue

        keys = ["DepMap_ID", "CCLE_Name", "COSMIC_ID", "Sanger ID"]
        for k in keys:
            if line[k] and not is_blank(line[k]):
                lookup[str(line[k])] = line["DepMap_ID"]

        if "MATCHED_NORMAL" in line["CCLE_Name"]:
            continue

        aliases = []
        prefix = line["CCLE_Name"].split('_')[0]
        blacklisted_prefixes = ["KMH2", "MS1", "TT"]
        if prefix not in blacklisted_prefixes:
            aliases.append(prefix)

        split_aliases = line["Aliases"].split(";")
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
