import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample, Aliquot, Individual, Project
from bmeg.edge import AliquotFor, BiosampleFor, InProject, PhenotypeOf
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from pydash import is_blank


def transform(path="source/ccle/CCLE_sample_info_file_2012-10-18.txt",
              prefix="ccle"):
    emitter = JSONEmitter(prefix)
    reader = bmeg.ioutils.read_tsv(path)

    # Load sample metadata.
    individual_gids = []
    project_gids = []
    phenotype_gids = []
    for row in reader:
        sample_id = row["CCLE name"]
        b = Biosample(sample_id, ccle_attributes=row)
        emitter.emit_vertex(b)

        a = Aliquot(aliquot_id=sample_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid(),
        )

        i = Individual(individual_id='CCLE:{}'.format(row["Cell line primary name"]),
                       ccle_attributes={'gender': row.get('Gender', None)})
        if i.gid() not in individual_gids:
            emitter.emit_vertex(i)
            individual_gids.append(i.gid())
        emitter.emit_edge(
            BiosampleFor(),
            b.gid(),
            i.gid(),
        )

        p = Project(project_id='CCLE:{}'.format(row["Site Primary"]))
        if p.gid() not in project_gids:
            emitter.emit_vertex(p)
            project_gids.append(p.gid())
        emitter.emit_edge(
            InProject(),
            i.gid(),
            p.gid(),
        )

        phenotype_name = row.get('Hist Subtype1', None)
        if not phenotype_name or is_blank(phenotype_name):
            phenotype_name = row.get('Histology')
        phenotype = phenotype_factory(name=phenotype_name)
        if phenotype.gid() not in phenotype_gids:
            emitter.emit_vertex(phenotype)
            phenotype_gids.append(phenotype.gid())
        emitter.emit_edge(
            PhenotypeOf(),
            a.gid(),
            phenotype.gid(),
        )

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
