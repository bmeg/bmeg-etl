import argparse
import os
import pandas as pd

from bmeg.vertex import COCACluster, Individual
from bmeg.edge import COCAClusterFor
from bmeg.emitter import JSONEmitter
from bmeg.utils import get_tcga_individual_barcode


def transform(args):
    """
    Hoadley et al. (2014). Cell-of-Origin Patterns Dominate the Molecular
    Classification of 10,000 Tumors from 33 Types of Cancer
    https://doi.org/10.1016/j.cell.2018.03.022

    Transform Table S1 downloaded from:
        https://ars.els-cdn.com/content/image/1-s2.0-S0092867418303027-mmc6.xlsx
    """
    emitter = JSONEmitter(args.output_prefix)

    df = pd.read_excel(args.input, sheet_name="Table S6 - iCluster", header=1)

    coca_clusters = {}
    for index, row in df.iterrows():
        individual_id = get_tcga_individual_barcode(row["Sample ID"])
        cluster_id = row["iCluster"]
        if cluster_id not in coca_clusters:
            # TODO add edge from each coca cluster to the pubmed vertex that
            # references this paper
            coca_clusters[cluster_id] = True
            coca = COCACluster(cluster_id=cluster_id)
            emitter.emit_vertex(coca)
        emitter.emit_edge(COCAClusterFor(),
                          from_gid=COCACluster.make_gid(cluster_id),
                          to_gid=Individual.make_gid(individual_id))

    emitter.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', "-i", type=str, required=True,
                        help='Path to the excel spreadsheet containing \
                        cluster assignments')
    parser.add_argument('--output-prefix', "-o", type=str, required=True,
                        help='Output file prefix')
    args = parser.parse_args()
    args.output_prefix = os.path.abspath(args.output_prefix)
    transform(args)
