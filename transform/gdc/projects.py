"""
Bulk download projects from GDC.
https://gdc.cancer.gov/
"""

from bmeg.util.cli import default_argument_parser
from bmeg.emitter import JSONEmitter
from bmeg import Project, Program

from transform.gdc.gdcutils import query_gdc, extract


parser = default_argument_parser()

expand_project_fields = ",".join("""
dbgap_accession_number
disease_type
name
primary_site
project_id
released
state
program
summary
""".strip().split())

keep_project_fields = """
dbgap_accession_number
disease_type
name
primary_site
summary
""".strip().split()


def transform(emitter):
    dedup = {}
    for row in query_gdc("projects", {"expand": expand_project_fields}):
        proj = Project(id=row["id"], **extract(row, keep_project_fields))
        prog = Program(id=row["program"]["name"], **row["program"])
        emitter.emit_vertex(proj)
        """
        emitter.emit_edge(
            HasProject(),
            to_gid=proj.gid(),
            from_gid=prog.gid(),
        )
        """
        if prog.id not in dedup:
            emitter.emit_vertex(prog)
            dedup[prog.id] = None


if __name__ == "__main__":
    args = parser.parse_args()
    emitter = JSONEmitter(directory="gdc")
    transform(emitter)
    emitter.close()
