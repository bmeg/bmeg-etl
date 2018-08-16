from bmeg.vertex import GeneExpression, Aliquot, ExpressionMetric
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id


p = "source/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct"
emitter = JSONEmitter("ccle")

for sample, values in parse_gct(p, "outputs/ccle", split_ensembl_id):
    g = GeneExpression(
        id=sample,
        source="ccle",
        scale=ExpressionMetric.RPKM,
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
