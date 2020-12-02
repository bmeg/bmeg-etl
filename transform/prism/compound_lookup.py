import json
import pandas

from bmeg.enrichers.drug_enricher import normalize


def create_table(primary_file="source/prism/primary-screen-replicate-collapsed-treatment-info.csv",
                 secondary_file="source/prism/secondary-screen-replicate-collapsed-treatment-info.csv"):
    df1 = pandas.read_csv(primary_file)
    df2 = pandas.read_csv(secondary_file)
    for drug in set(df1.name.tolist() + df2.name.tolist()):
        if pandas.isna(drug):
            print(drug, "", sep="\t")
            continue
        resp = normalize(drug)
        if resp:
            resp = json.dumps(resp)
        else:
            resp = ""
        print(drug, resp, sep="\t")


if __name__ == "__main__":
    create_table()
