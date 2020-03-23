import json
import pandas

from bmeg.enrichers.drug_enricher import normalize


def create_table(in_file="source/prism/primary-screen-replicate-collapsed-treatment-info.csv"):
    df = pandas.read_csv(in_file, index_col=0)
    for drug in set(df.name):
        resp = normalize(drug)
        if resp:
            resp = json.dumps(resp)
        else:
            resp = ""
        print(drug, resp, sep="\t")


if __name__ == "__main__":
    create_table()
