from bmeg.vertex import Expression, Aliquot, ExpressionMetric
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id
from bmeg.utils import ensure_directory


def transform(path="source/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct",
              emitter_prefix=None,
              emitter_directory="ccle"):
    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    # parse_gct creates a large file 'tmp_matrix'  which we don't want to put in /tmp
    # at the same time, we need to ensure an output path exists (ie. for travis)
    ensure_directory("outputs/ccle")
    for sample, values in parse_gct(path, "outputs/ccle", split_ensembl_id):
        # strip out broad suffix
        sample = sample.split()[0]
        g = Expression(
            id=sample,
            source="ccle",
            metric=ExpressionMetric.RPKM,
            method="Unknown",
            values=values,
        )
        emitter.emit_vertex(g)
        emitter.emit_edge(
            ExpressionOf(),
            from_gid=g.gid(),
            to_gid=Aliquot.make_gid(sample)
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
