import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample, Aliquot, Individual, Project
from bmeg.edge import AliquotFor, BiosampleFor, InProject, PhenotypeOf
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from pydash import is_blank


def transform(path="source/ccle/DepMap-2018q3-celllines.csv",
              prefix="ccle"):
    emitter = JSONEmitter(prefix)
    reader = bmeg.ioutils.read_csv(path)

    # Load sample metadata.
    individual_gids = []
    project_gids = []
    phenotype_gids = []
    # Broad_ID,CCLE_Name,Aliases,COSMIC_ID,Sanger ID,Primary Disease,Subtype Disease,Gender,Source
    # ACH-000557,AML193_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE,AML-193,,,Leukemia,,Female,ATCC
    for row in reader:
        sample_id = row["Broad_ID"]
        b = Biosample(sample_id, ccle_attributes=row)
        emitter.emit_vertex(b)

        a = Aliquot(aliquot_id=sample_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid(),
        )

        i = Individual(individual_id='CCLE:{}'.format(sample_id),
                       ccle_attributes={'gender': row.get('Gender', None)})
        if i.gid() not in individual_gids:
            emitter.emit_vertex(i)
            individual_gids.append(i.gid())
        emitter.emit_edge(
            BiosampleFor(),
            b.gid(),
            i.gid(),
        )

        project_id = row["CCLE_Name"]
        p = Project(project_id='CCLE:{}'.format(project_id))
        if p.gid() not in project_gids:
            emitter.emit_vertex(p)
            project_gids.append(p.gid())
        emitter.emit_edge(
            InProject(),
            i.gid(),
            p.gid(),
        )

        phenotype_name = row.get('Subtype Disease', None)
        if not phenotype_name or is_blank(phenotype_name):
            phenotype_name = row.get('Primary Disease')
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
