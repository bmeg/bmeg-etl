from bmeg.vertex import Expression, Aliquot, ExpressionMetric
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id


p = "source/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"
emitter = JSONEmitter("gtex")

for sample, values in parse_gct(p, "outputs/gtex", split_ensembl_id):
    g = Expression(
        id=sample,
        source="gtex",
        metric=ExpressionMetric.GENE_TPM,
        method="Illumina HiSeq",
        values=values,
    )
    emitter.emit_vertex(g)
    emitter.emit_edge(
        ExpressionOf(),
        from_gid=g.gid(),
        to_gid=Aliquot.make_gid(sample)
    )

emitter.close()
