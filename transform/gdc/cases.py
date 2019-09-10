"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""
import json
import pydash

from bmeg.util.cli import default_argument_parser
from bmeg.emitter import JSONEmitter
from bmeg import (Sample, Aliquot, Case, Project, Program,
                  Aliquot_Sample_Sample,
                  Sample_Case_Case,
                  Case_Projects_Project,
                  Sample_Projects_Project,
                  Aliquot_Projects_Project,
                  Project_Programs_Program,
                  Case_Phenotypes_Phenotype,
                  Sample_Phenotypes_Phenotype)
from bmeg.enrichers.phenotype_enricher import phenotype_factory
from transform.gdc.gdcutils import extract


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

    phenotypes = {}
    programs = {}
    projects = {}
    with open(input_path, "r") as fh:
        for row in fh:
            row = json.loads(row)

            # only emit TCGA data
            if pydash.get(row, "project.program.name", "") != "TCGA":
                continue

            # program
            prog = Program(id=Program.make_gid(row["project"]["program"]["name"]),
                           program_id=row["project"]["program"]["name"],
                           gdc_attributes=extract(row["project"]["program"], keep_program_fields))
            if prog.gid() not in programs:
                emitter.emit_vertex(prog)
                programs[prog.gid()] = True

            # project
            proj = Project(id=Project.make_gid(row["project"]["project_id"]),
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
            c = Case(id=Case.make_gid(row["id"]),
                     submitter_id=row["submitter_id"],
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

            # create phenotype
            pheno_name = pydash.get(row, "project.project_id", "unknown").replace("TCGA-", "")
            pheno = phenotype_factory(pheno_name)
            if pheno.gid() not in phenotypes:
                emitter.emit_vertex(pheno)
                phenotypes[pheno.gid()] = None
            # case <-> phenotype edges
            emitter.emit_edge(
                Case_Phenotypes_Phenotype(
                    from_gid=c.gid(),
                    to_gid=pheno.gid()
                ),
                emit_backref=True
            )

            for sample in row.get("samples", []):
                sample_fields = extract(
                    sample,
                    ["tumor_descriptor", "sample_type", "submitter_id"],
                )
                # sample
                s = Sample(id=Sample.make_gid(sample["sample_id"]),
                           submitter_id=sample["submitter_id"],
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
                if "Normal" not in sample_fields["sample_type"]:
                    # sample <-> phenotype edges
                    emitter.emit_edge(
                        Sample_Phenotypes_Phenotype(
                            from_gid=s.gid(),
                            to_gid=pheno.gid()
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
                            a = Aliquot(id=Aliquot.make_gid(aliquot["aliquot_id"]),
                                        submitter_id=aliquot["submitter_id"],
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
