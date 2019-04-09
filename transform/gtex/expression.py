from bmeg.vertex import GeneExpression, Aliquot, ExpressionMetric
from bmeg.edge import HasGeneExpression
from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id


p = "source/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
emitter = JSONEmitter("gtex")

for sample, values in parse_gct(p, "outputs/gtex", split_ensembl_id):
    g = GeneExpression(
        id=sample,
        source="gtex",
        metric=ExpressionMetric.GENE_TPM,
        method="Illumina HiSeq",
        values=values,
    )
    emitter.emit_vertex(g)
    emitter.emit_edge(
        HasGeneExpression(),
        to_gid=g.gid(),
        from_gid=Aliquot.make_gid(sample)
    )

emitter.close()
