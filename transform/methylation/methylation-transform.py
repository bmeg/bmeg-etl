import argparse
import os

import bmeg.ioutils
from bmeg.models.vertex_models import Callset, Gene, MethylationProbe
from bmeg.models.edge_models import (MethlyationProbeValue,
                                     MethlyationProbeFor)
from bmeg.models.conversions import query_mygeneinfo_for_ensembl_gene
from bmeg.models.emitter import Emitter
from bmeg.models.utils import get_tcga_sample_barcode, tcga_barcode_is_tumor


def transform(args):
    """
    Transform legacy methylation data from the gdc. Example:
        https://portal.gdc.cancer.gov/legacy-archive/files/85d8d46e-5206-4c07-95e0-56b22c013321
    """
    emitter = Emitter(args.output_prefix)
    emit = emitter.emit

    reader = bmeg.ioutils.tsv(args.input)

    sample_id = get_tcga_sample_barcode(reader.fieldsnames[1])
    if not tcga_barcode_is_tumor(sample_id):
        raise RuntimeError("%s is not a tumor barcode" % sample_id)

    for line in reader:
        gene = line["Gene_Symbol"]
        gene = query_mygeneinfo_for_ensembl_gene(gene)

        p = MethylationProbe(genome=args.genome,
                             chromosome=line["Chromosome"],
                             start=int(line["Genomic_Coordinate"]),
                             end=int(line["Genomic_Coordinate"]),
                             probe_id=line["Composite Element REF"])

        c = Callset(tumor_biosample_id=sample_id,
                    normal_biosample_id="",
                    call_method="Illumina Human Methylation 450",
                    source="GDC")

        mpv = MethlyationProbeValue(from_gid=p.gid,
                                    to_gid=c.gid,
                                    value=line["Beta_value"])

        mpfg = MethlyationProbeFor(from_gid=p.gid,
                                   to_gid=Gene.make_gid(gene))

        emit(p)
        emit(c)
        emit(mpv)
        emit(mpfg)

    emitter.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', "-i", type=str, required=True,
                        help='Path to DNA methylation text file')
    parser.add_argument('--output-prefix', "-o", type=str, required=True,
                        help='Output file prefix')
    parser.add_argument('--genome', "-g", type=str, default="GRCh37",
                        help='Reference genome build')
    args = parser.parse_args()
    args.output_prefix = os.path.abspath(args.output_prefix)
    transform(args)
