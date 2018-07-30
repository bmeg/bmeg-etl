import argparse
import os

import bmeg.ioutils
from bmeg.vertex import Callset, Gene
from bmeg.edge import GeneCNAValueCall
from bmeg.conversions import ensembl_gene_lookup
from bmeg.emitter import JSONEmitter
from bmeg.utils import get_tcga_sample_barcode, tcga_barcode_is_tumor


def transform(args):
    """
    Transform the file downloaded from:
        https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
    """
    emitter = JSONEmitter(args.output_prefix)
    emit = emitter.emit

    for line in bmeg.ioutils.tsv(args.input):
        identifier = line["Sample"]
        del line["Sample"]

        ids = identifier.split("|")
        if len(ids) == 2:
            gene = ids[1].split(".")[0]
        else:
            gene = ids[0]
            gene = ensembl_gene_lookup(gene)

        for id, val in line.items():
            sample_id = get_tcga_sample_barcode(id)
            if not tcga_barcode_is_tumor(sample_id):
                raise RuntimeError("%s is not a tumor barcode" % sample_id)
            c = Callset(tumor_biosample_id=sample_id,
                        normal_biosample_id="",
                        call_method="GISTIC2",
                        source="UCSCXena")
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
    parser.add_argument('--output-prefix', "-o", type=str, required=True,
                        help='Output file prefix')
    args = parser.parse_args()
    args.output_prefix = os.path.abspath(args.output_prefix)
    transform(args)
