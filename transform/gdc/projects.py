"""
Bulk download projects from GDC.
https://gdc.cancer.gov/
"""

from bmeg.util.cli import default_argument_parser
from bmeg.edge import InProject, BiosampleFor, AliquotFor
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Project

from gdcutils import query_gdc, extract


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
disease_type
name
primary_site
program
summary
""".strip().split()


def transform(emitter):
    for row in query_gdc("projects", {"expand": expand_project_fields}):
        p = Project(row["id"], extract(row, keep_project_fields))
        emitter.emit_vertex(p)


if __name__ == "__main__":
    args = parser.parse_args()
    emitter = JSONEmitter("gdc")
    transform(emitter)
    emitter.close()
