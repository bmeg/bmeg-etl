from bmeg.ioutils import read_tsv


def extract_case_id(sample_id):
    return "-".join(sample_id.split("-")[:2])


def create_lookups(samples_path="source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt",
                   output_dir="source/gtex"):

    samples = read_tsv(samples_path)

    project_lookup = {}
    for row in samples:
        sample_id = row["SAMPID"]
        case_id = extract_case_id(sample_id)
        proj_id = "GTEx_{}".format(row["SMTS"])

        project_lookup[sample_id] = proj_id
        project_lookup[case_id] = proj_id

    with open("{}/project_lookup.tsv".format(output_dir), "w") as handle:
        for key, value in project_lookup.items():
            handle.write("{}\t{}\n".format(key, value))

    return


if __name__ == "__main__":
    create_lookups()
