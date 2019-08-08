import json


def create_lookups(input_path="source/gdc/cases.json",
                   output_dir="source/gdc"):

    id_lookup = {}
    project_lookup = {}
    with open(input_path, "r") as fh:
        for row in fh:
            row = json.loads(row)

            # project
            proj_id = row["project"]["project_id"]

            # case
            project_lookup[row["id"]] = proj_id
            project_lookup[row["submitter_id"]] = proj_id

            id_lookup[row["id"]] = row["id"]
            id_lookup[row["submitter_id"]] = row["id"]

            for sample in row.get("samples", []):
                # sample
                project_lookup[sample["sample_id"]] = proj_id
                project_lookup[sample["submitter_id"]] = proj_id

                id_lookup[sample["sample_id"]] = sample["sample_id"]
                id_lookup[sample["submitter_id"]] = sample["sample_id"]

                for portion in sample.get("portions", []):
                    for analyte in portion.get("analytes", []):
                        for aliquot in analyte.get("aliquots", []):
                            project_lookup[aliquot["aliquot_id"]] = proj_id
                            project_lookup[aliquot["submitter_id"]] = proj_id

                            id_lookup[aliquot["aliquot_id"]] = aliquot["aliquot_id"]
                            id_lookup[aliquot["submitter_id"]] = aliquot["aliquot_id"]

    with open("{}/project_lookup.tsv".format(output_dir), "w") as handle:
        for key, value in project_lookup.items():
            handle.write("{}\t{}\n".format(key, value))

    with open("{}/id_lookup.tsv".format(output_dir), "w") as handle:
        for key, value in id_lookup.items():
            handle.write("{}\t{}\n".format(key, value))
    return


if __name__ == "__main__":
    create_lookups()
