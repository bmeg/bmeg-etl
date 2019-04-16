import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg import Sample, Aliquot, Case, Project, Program
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from pydash import is_blank


def transform(path="source/ccle/DepMap-2018q4-celllines.csv",
              prefix="ccle"):
    emitter = JSONEmitter(directory=prefix)
    reader = bmeg.ioutils.read_csv(path)

    # Load sample metadata.
    case_gids = []
    project_gids = []
    phenotype_gids = []
    prog = Program(id="DepMap", name="DepMap", dbgap_accession_number="NA")
    emitter.emit_vertex(prog)
    # DepMap_ID,CCLE_Name,Aliases,COSMIC_ID,Sanger ID,Primary Disease,Subtype Disease,Gender,Source
    # ACH-000001,NIHOVCAR3_OVARY,NIH:OVCAR-3;OVCAR3,905933,2201,Ovarian Cancer,"Adenocarcinoma, high grade serous",Female,ATCC
    for row in reader:
        sample_id = row["DepMap_ID"]
        s = Sample(id=sample_id, submitter_id=row["CCLE_Name"], **row)
        emitter.emit_vertex(s)

        a = Aliquot(id=sample_id)
        emitter.emit_vertex(a)
        #emitter.emit_edge(
        #    HasAliquot(),
        #    to_gid=a.gid(),
        #    from_gid=s.gid()
        #)

        i = Case(id=sample_id, gender=row.get('Gender', None))
        if i.id not in case_gids:
            emitter.emit_vertex(i)
            case_gids.append(i.id)
        #emitter.emit_edge(
        #    HasSample(),
        #    to_gid=s.gid(),
        #    from_gid=i.gid()
        #)

        project_id = "DepMap_{}".format("_".join(row["Primary Disease"].split()))
        # create project
        proj = Project(id=project_id)
        if proj.id not in project_gids:
            emitter.emit_vertex(proj)
            #emitter.emit_edge(
            #    HasProject(),
            #    to_gid=proj.gid(),
            #    from_gid=prog.gid()
            #)
            project_gids.append(proj.id)
        #emitter.emit_edge(
        #    HasCase(),
        #    to_gid=i.gid(),
        #    from_gid=proj.gid()
        #)

        phenotype_name = row.get('Subtype Disease', None)
        if not phenotype_name or is_blank(phenotype_name):
            phenotype_name = row.get('Primary Disease')
        phenotype = phenotype_factory(name=phenotype_name)
        if phenotype.term_id not in phenotype_gids:
            emitter.emit_vertex(phenotype)
            phenotype_gids.append(phenotype.term_id)
        #emitter.emit_edge(
        #    PhenotypeOf(),
        #    a.gid(),
        #    phenotype.gid(),
        #)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
