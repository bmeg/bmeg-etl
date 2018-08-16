from bmeg.vertex import GeneExpression, Aliquot, ExpressionMetric
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id


p = "source/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"
emitter = JSONEmitter("gtex")

for sample, values in parse_gct(p, "outputs/gtex", split_ensembl_id):
    g = GeneExpression(
        id=sample,
        source="gtex",
        scale=ExpressionMetric.TPM,
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
