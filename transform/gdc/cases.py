"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""
import json

from bmeg.util.cli import default_argument_parser
from bmeg.emitter import JSONEmitter
from bmeg import (Sample, Aliquot, Case, Project, Program,
                  Aliquot_Sample_Sample,
                  Sample_Case_Case,
                  Case_Projects_Project,
                  Sample_Projects_Project,
                  Aliquot_Projects_Project,
                  Project_Programs_Program)
from gdcutils import extract


# These are the fields we want to keep from the GDC Case (BMEG Case).
keep_case_fields = """
diagnoses
demographic
disease_type
primary_site
summary
""".strip().split()

keep_project_fields = """
dbgap_accession_number
disease_type
name
primary_site
summary
""".strip().split()

keep_program_fields = """
dbgap_accession_number
name
""".strip().split()


def transform(input_path="source/gdc/cases.json",
              emitter_prefix=None,
              emitter_directory="gdc"):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    programs = {}
    projects = {}
    with open(input_path, "r") as fh:
        for row in fh:
            row = json.loads(row)

            # program
            prog = Program(submitter_id=Program.make_gid(row["project"]["program"]["name"]),
                           program_id=row["project"]["program"]["name"],
                           gdc_attributes=extract(row["project"]["program"], keep_program_fields))
            if prog.gid() not in programs:
                emitter.emit_vertex(prog)
                programs[prog.gid()] = True

            # project
            proj = Project(submitter_id=Project.make_gid(row["project"]["project_id"]),
                           project_id=row["project"]["project_id"],
                           gdc_attributes=extract(row["project"], keep_project_fields))
            if proj.gid() not in projects:
                emitter.emit_vertex(proj)
                projects[proj.gid()] = True
                # project <-> program edges
                emitter.emit_edge(
                    Project_Programs_Program(
                        from_gid=proj.gid(),
                        to_gid=prog.gid()
                    ),
                    emit_backref=True
                )
            project_gid = proj.gid()

            # case
            c = Case(submitter_id=Case.make_gid(row["id"]),
                     case_id=row["id"],
                     gdc_attributes=extract(row, keep_case_fields),
                     project_id=project_gid)
            emitter.emit_vertex(c)
            # case <-> project edges
            emitter.emit_edge(
                Case_Projects_Project(
                    from_gid=c.gid(),
                    to_gid=project_gid,
                ),
                emit_backref=True
            )

        for sample in row.get("samples", []):
            sample_fields = extract(
                sample,
                ["tumor_descriptor", "sample_type", "submitter_id"],
            )
            # sample
            s = Sample(submitter_id=Sample.make_gid(sample["sample_id"]),
                       sample_id=sample["sample_id"],
                       gdc_attributes=sample_fields,
                       project_id=project_gid)
            emitter.emit_vertex(s)
            # sample <-> case edges
            emitter.emit_edge(
                Sample_Case_Case(
                    from_gid=s.gid(),
                    to_gid=c.gid()
                ),
                emit_backref=True
            )
            # sample <-> project edges
            emitter.emit_edge(
                Sample_Projects_Project(
                    from_gid=s.gid(),
                    to_gid=project_gid
                ),
                emit_backref=True
            )

            for portion in sample.get("portions", []):
                for analyte in portion.get("analytes", []):
                    for aliquot in analyte.get("aliquots", []):
                        aliquot_fields = extract(
                            aliquot,
                            ["analyte_type", "submitter_id", "aliquot_id"],
                        )
                        fields = dict(sample_fields)
                        fields.update(aliquot_fields)
                        # aliquot
                        a = Aliquot(submitter_id=Aliquot.make_gid(aliquot["aliquot_id"]),
                                    aliquot_id=aliquot["aliquot_id"],
                                    gdc_attributes=fields,
                                    project_id=project_gid)
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
                                to_gid=project_gid
                            ),
                            emit_backref=True
                        )
    emitter.close()


if __name__ == "__main__":
    parser = default_argument_parser()
    args = parser.parse_args()
    transform()
