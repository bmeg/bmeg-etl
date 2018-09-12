from bmeg.vertex import Expression, Aliquot, ExpressionMetric
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id


def transform():

    p = "source/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct"
    emitter = JSONEmitter("ccle")

    for sample, values in parse_gct(p, "outputs/ccle", split_ensembl_id):
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
