import argparse
import os
import pandas as pd

from bmeg.models.vertex_models import COCACluster, Individual
from bmeg.models.edge_models import COCAClusterContains
from bmeg.models.emitter import Emitter
from bmeg.models.utils import get_tcga_individual_barcode


def transform(args):
    """
    Hoadley et al. (2014). Multiplatform Analysis of 12 Cancer Types Reveals
    Molecular Classification within and across Tissues of Origin
    https://doi.org/10.1016/j.cell.2014.06.049

    Transform Table S1 downloaded from:
        https://ars.els-cdn.com/content/image/1-s2.0-S0092867414008769-mmc2.xlsx
    """
    emitter = Emitter(args.output_prefix)
    emit = emitter.emit

    df = pd.read_excel(args.input, sheet_name="Sheet1", header=1)

    coca_clusters = {}
    for index, row in df.iterrows():
        individual_id = get_tcga_individual_barcode(row["Sample"])
        cluster_id = row["COCA"]
        if cluster_id not in coca_clusters:
            coca_clusters[cluster_id] = True
            coca = COCACluster(cluster_id=cluster_id)
            emit(coca)
        cvfg = COCAClusterContains(from_gid=COCACluster.make_gid(cluster_id),
                                   to_gid=Individual.make_gid(individual_id))
        emit(cvfg)
    emitter.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', "-i", type=str, required=True,
                        help='Path to the excel spreadsheet containing COCA \
                        cluster assignments')
    parser.add_argument('--output-prefix', "-o", type=str, required=True,
                        help='Output file prefix')
    args = parser.parse_args()
    args.output_prefix = os.path.abspath(args.output_prefix)
    transform(args)
