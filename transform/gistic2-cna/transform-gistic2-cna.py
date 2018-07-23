import argparse
import csv
import gzip
import io
import os

from bmeg.models.vertex_models import Callset, Gene
from bmeg.models.edge_models import GeneCNAValueCall
from bmeg.models.conversions import query_mygeneinfo_for_ensembl_gene
from bmeg.models.emitter import Emitter
from bmeg.models.utils import get_tcga_sample_barcode, tcga_barcode_is_tumor


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
        identifier = line["Sample"]
        del line["Sample"]

        ids = identifier.split("|")
        if len(ids) == 2:
            gene = ids[1].split(".")[0]
        else:
            gene = ids[0]
            gene = query_mygeneinfo_for_ensembl_gene(gene)

        for id, val in line.items():
            sample_id = get_tcga_sample_barcode(id)
            if not tcga_barcode_is_tumor(sample_id):
                raise RuntimeError("%s is not a tumor barcode" % sample_id)
            c = Callset(tumor_biosample_id=sample_id,
                        normal_biosample_id="",
                        call_method="GISTIC2")
            cvfg = GeneCNAValueCall(from_gid=c.gid,
                                    to_gid=Gene.make_gid(gene),
                                    value=val)
            emit(c)
            emit(cvfg)
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
