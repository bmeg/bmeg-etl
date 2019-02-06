import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample, Aliquot, Individual, Project
from bmeg.edge import AliquotFor, BiosampleFor, InProject, PhenotypeOf
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from pydash import is_blank

PROJECT_CORRECTIONS = {
    # suffix changes
    'TESTES': 'TESTIS',
    'HAEMATOPOIETIC_AND_LYMPHOID': 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE',
    # wholesale name changes
    '[MERGED_TO_ACH-000474]NCIH292_LUNG': 'LUNG',
    '[MERGED_TO_ACH-000109]H3255_LUNG': 'LUNG'
}


def transform(path="source/ccle/DepMap-2018q4-celllines.csv",
              prefix="ccle"):
    emitter = JSONEmitter(directory=prefix)
    reader = bmeg.ioutils.read_csv(path)

    # Load sample metadata.
    individual_gids = []
    project_gids = []
    phenotype_gids = []
    # DepMap_ID,CCLE_Name,Aliases,COSMIC_ID,Sanger ID,Primary Disease,Subtype Disease,Gender,Source
    # ACH-000001,NIHOVCAR3_OVARY,NIH:OVCAR-3;OVCAR3,905933,2201,Ovarian Cancer,"Adenocarcinoma, high grade serous",Female,ATCC
    for row in reader:
        sample_id = row["DepMap_ID"]
        sample_id = 'CCLE:{}'.format(sample_id)
        b = Biosample(biosample_id=sample_id,
                      ccle_attributes=row)
        emitter.emit_vertex(b)

        a = Aliquot(aliquot_id=sample_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid(),
        )

        i = Individual(individual_id=sample_id,
                       ccle_attributes={'gender': row.get('Gender', None)})
        if i.gid() not in individual_gids:
            emitter.emit_vertex(i)
            individual_gids.append(i.gid())
        emitter.emit_edge(
            BiosampleFor(),
            b.gid(),
            i.gid(),
        )

        # first see if we have wholesale name changes
        project_id = PROJECT_CORRECTIONS.get(row["CCLE_Name"], row["CCLE_Name"])
        # strip off prefix
        name_parts = project_id.split('_')
        name_start = 1
        if len(name_parts) == 1:
            name_start = 0
        project_id = '_'.join(name_parts[name_start:])
        # suffix name changes
        project_id = PROJECT_CORRECTIONS.get(project_id, project_id)
        # create project
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
