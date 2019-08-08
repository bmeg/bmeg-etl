from bmeg.emitter import JSONEmitter
from bmeg.ioutils import read_tsv
from bmeg import (Sample, Aliquot, Case, Project, Program,
                  Aliquot_Sample_Sample,
                  Sample_Case_Case,
                  Case_Projects_Project,
                  Sample_Projects_Project,
                  Aliquot_Projects_Project,
                  Project_Programs_Program)


def extract_case_id(sample_id):
    return "-".join(sample_id.split("-")[:2])


def transform(cases_path="source/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt",
              samples_path="source/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt",
              emitter_prefix="gtex",
              emitter_directory="gtex"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)
    samples = read_tsv(samples_path)
    cases = read_tsv(cases_path)

    prog = Program(id=Program.make_gid("GTEx"),
                   program_id="GTEx")
    emitter.emit_vertex(prog)

    case_props = {}
    for row in cases:
        case_id = row["SUBJID"]
        case_props[case_id] = row

    case_ids = {}
    project_ids = {}
    for row in samples:
        sample_id = row["SAMPID"]
        case_id = extract_case_id(sample_id)

        # project
        proj_name = "GTEx_{}".format(row["SMTS"])
        proj = Project(id=Project.make_gid(proj_name),
                       project_id=proj_name)
        if proj.gid() not in project_ids:
            emitter.emit_vertex(proj)
            # project <-> program edges
            emitter.emit_edge(
                Project_Programs_Program(
                    from_gid=proj.gid(),
                    to_gid=prog.gid()
                ),
                emit_backref=True
            )
            project_ids[proj.gid()] = True

        # case
        c = Case(
            id=Case.make_gid(case_id),
            submitter_id=case_id,
            case_id=case_id,
            gtex_attributes=case_props[case_id],
            project_id=proj.gid()
        )
        if c.gid() not in case_ids:
            emitter.emit_vertex(c)
            # case <-> project edges
            emitter.emit_edge(
                Case_Projects_Project(
                    from_gid=Case.make_gid(case_id),
                    to_gid=proj.gid(),
                ),
                emit_backref=True
            )
            case_ids[c.gid()] = True

        # sample
        s = Sample(
            id=Sample.make_gid(sample_id),
            submitter_id=sample_id,
            sample_id=sample_id,
            gtex_attributes=row,
            project_id=proj.gid()
        )
        emitter.emit_vertex(s)
        # sample <-> case edges
        emitter.emit_edge(
            Sample_Case_Case(
                from_gid=s.gid(),
                to_gid=Case.make_gid(case_id)
            ),
            emit_backref=True
        )
        # sample <-> project edges
        emitter.emit_edge(
            Sample_Projects_Project(
                from_gid=s.gid(),
                to_gid=proj.gid()
            ),
            emit_backref=True
        )

        # aliquot
        a = Aliquot(
            id=Aliquot.make_gid(sample_id),
            submitter_id=sample_id,
            aliquot_id=sample_id,
            project_id=proj.gid()
        )
        emitter.emit_vertex(a)
        # aliquot <-> sample edges
        emitter.emit_edge(
            Aliquot_Sample_Sample(
                from_gid=a.gid(),
                to_gid=s.gid()
            ),
            emit_backref=True
        )
        # aliquot <-> project edges
        emitter.emit_edge(
            Aliquot_Projects_Project(
                from_gid=a.gid(),
                to_gid=proj.gid()
            ),
            emit_backref=True
        )

    emitter.close()


if __name__ == "__main__":
    transform()
