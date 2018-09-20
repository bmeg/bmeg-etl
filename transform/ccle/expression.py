from bmeg.vertex import Expression, Aliquot, ExpressionMetric
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id


def transform(path="source/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct",
              emitter_prefix="ccle",
              emitter_directory="outputs/ccle"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

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
