"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""

import json

import bmeg.models.emitter
from bmeg.models.vertex_models import Individual, Biosample, Project

from gdcutils import query_gdc, extract


emitter = bmeg.models.emitter.JSONEmitter("gdc")

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

for row in query_gdc("projects", {"expand": expand_project_fields}):
    p = Project(row["id"], extract(row, keep_project_fields))
    emitter.emit(p)

emitter.close()
