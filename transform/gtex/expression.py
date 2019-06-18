from bmeg import (GeneExpression, Aliquot, Project,
                  GeneExpression_Aliquot_Aliquot)

from bmeg.emitter import JSONEmitter
from bmeg.gct import parse_gct, split_ensembl_id
from bmeg.ioutils import read_lookup


def transform(gct_file="source/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz",
              project_lookup_path="source/gtex/project_lookup.tsv",
              emitter_prefix=None,
              emitter_directory="gtex"):

    projects = read_lookup(project_lookup_path)

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    for sample, values in parse_gct(gct_file, "outputs/gtex", split_ensembl_id):
        g = GeneExpression(
            submitter_id=GeneExpression.make_gid(sample),
            metric="GENE_TPM",
            method="Illumina HiSeq",
            values=values,
            project_id=Project.make_gid(projects.get(sample, None))
        )
        emitter.emit_vertex(g)
        emitter.emit_edge(
            GeneExpression_Aliquot_Aliquot(
                from_gid=g.gid(),
                to_gid=Aliquot.make_gid(sample)
            ),
            emit_backref=True
        )

    emitter.close()


if __name__ == "__main__":
    transform()
