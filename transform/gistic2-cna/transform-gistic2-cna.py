import argparse
import csv
import gzip
import io
import os

from bmeg.models.vertex_models import Callset, Gene
from bmeg.models.edge_models import CNAValueForGene
from bmeg.models.conversions import query_mygeneinfo_for_ensembl_gene
from bmeg.models.emitter import Emitter


def transform(args):
    """
    Transform the file downloaded from:
        https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
    """
    emitter = Emitter(args.output_prefix)
    emit = emitter.emit

    if args.gz:
        inhandle = io.TextIOWrapper(gzip.GzipFile(args.input))
    else:
        inhandle = open(args.input, "r")

    reader = csv.DictReader(inhandle, delimiter="\t")
    for line in reader:
        gene = line["Sample"]
        del line["Sample"]
        ensembl_gene_id = query_mygeneinfo_for_ensembl_gene(gene)
        for k, val in line.items():
            vert = Callset(tumor_biosample_id=k,
                           normal_biosample_id="",
                           call_method="GISTIC2")
            edge = CNAValueForGene(from_gid=Gene.make_gid(ensembl_gene_id),
                                   to_gid=vert.gid,
                                   value=val)
            emit(vert)
            emit(edge)
        emitter.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', "-i", type=str, required=True,
                        help='Path to the GISTIC2 copy number matrix file')
    parser.add_argument('--gz', action="store_true", default=False,
                        help='Input file is gzipped')
    parser.add_argument('--output-prefix', "-o", type=str, required=True,
                        help='Output file prefix')
    args = parser.parse_args()
    args.output_prefix = os.path.abspath(args.output_prefix)
    transform(args)
