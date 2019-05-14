from bmeg.vertex import GeneExpression, TranscriptExpression, Aliquot, ExpressionMetric
from bmeg.edge import GeneExpressionOf, TranscriptExpressionOf
from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id
from bmeg.utils import ensure_directory

import os
import re
import pandas


def transform_rpkm(path="source/depmap/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct",
                   emitter_prefix=None,
                   emitter_directory="depmap"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    outdir = os.path.join("outputs", emitter_directory)
    # parse_gct creates a large file 'tmp_matrix'  which we don't want to put in /tmp
    # at the same time, we need to ensure an output path exists (ie. for travis)
    ensure_directory(outdir)
    for sample, values in parse_gct(path, outdir, split_ensembl_id):
        # use broad suffix  "QGP1_PANCREAS (ACH-000347)" -> "ACH-000347"
        sample = sample.split()[1].replace('(', '').replace(')', '')
        g = GeneExpression(
            id=sample,
            source="ccle",
            metric=ExpressionMetric.RPKM,
            method="Unknown",
            values=values,
        )
        emitter.emit_vertex(g)
        emitter.emit_edge(
            GeneExpressionOf(),
            from_gid=g.gid(),
            to_gid=Aliquot.make_gid(sample)
        )

    emitter.close()


def transform_gene_tpm(path="source/ccle/CCLE_depMap_19Q1_TPM.csv",
                       emitter_prefix=None,
                       emitter_directory="ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    outdir = os.path.join("outputs", emitter_directory)
    ensure_directory(outdir)

    tmp_df = pandas.read_csv(path, sep=",", index_col=0)
    tmp_df.rename(columns=lambda x: re.search(r'\((ENSG.*)\)', x).group(1), inplace=True)
    for sample, values in tmp_df.iterrows():
        g = GeneExpression(
            id=sample,
            source="CCLE",
            metric=ExpressionMetric.GENE_TPM,
            method="Unknown",
            values=values.to_dict(),
        )
        emitter.emit_vertex(g)
        emitter.emit_edge(
            GeneExpressionOf(),
            from_gid=g.gid(),
            to_gid=Aliquot.make_gid("CCLE:%s:GeneExpression" % (sample))
        )

    emitter.close()


def transform_tpm(path="source/ccle/CCLE_depMap_19Q1_TPM_transcripts.csv",
                  emitter_prefix=None,
                  emitter_directory="ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    outdir = os.path.join("outputs", emitter_directory)
    ensure_directory(outdir)

    tmp_df = pandas.read_csv(path, sep=",", index_col=0)
    tmp_df.rename(columns=lambda x: re.search(r'\((ENST.*)\)', x).group(1), inplace=True)
    for sample, values in tmp_df.iterrows():
        t = TranscriptExpression(
            id=sample,
            source="CCLE",
            metric=ExpressionMetric.TPM,
            method="Unknown",
            values=values.to_dict(),
        )
        emitter.emit_vertex(t)
        emitter.emit_edge(
            TranscriptExpressionOf(),
            from_gid=t.gid(),
            to_gid=Aliquot.make_gid("CCLE:%s:TranscriptExpression" % (sample))
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform_tpm()
    transform_gene_tpm()
