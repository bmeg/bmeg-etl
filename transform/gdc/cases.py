"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""

import json

import bmeg.models.emitter
from bmeg.models.vertex_models import Individual, Biosample, Project

from gdcutils import extract, query_gdc

#emitter = JSONEmitter("gdc")
#emitter = BSONEmitter("gdc")
#emitter = MsgpackEmitter("gdc")
emitter = bmeg.models.emitter.DebugEmitter()

# The GDC API requires you to explicitly request that nested fields be expanded.
# https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/#cases-field-groups
#
# Note that (as of this writing) we are expanding most but not all possible fields.
# Mostly we're skipping "files" data.
expand_case_fields = ",".join("""
demographic
diagnoses
diagnoses.treatments
exposures
family_histories
project
project.program
samples
samples.annotations
samples.portions
samples.portions.analytes
samples.portions.analytes.aliquots
samples.portions.analytes.aliquots.annotations
samples.portions.analytes.aliquots.center
samples.portions.analytes.annotations
samples.portions.annotations
samples.portions.center
samples.portions.slides
samples.portions.slides.annotations
summary
summary.data_categories
summary.experimental_strategies
tissue_source_site
""".strip().split())

# These are the fields we want to keep from the GDC Case (BMEG Individual).
keep_case_fields = """
diagnoses
demographic
disease_type
primary_site
summary
project
""".strip().split()

# Crawl all cases, samples, aliquots to generate BMEG Individuals and Biosamples.
for row in query_gdc("cases", {"expand": expand_case_fields}):

    i = Individual(row["id"], extract(row, keep_case_fields))
    emitter.emit(i)

    for sample in row.get("samples", []):
        print(json.dumps(sample, indent=True))

        sample_fields = extract(sample, ["tumor_descriptor", "sample_type", "submitter_id"])
        b = Biosample(sample["sample_id"], sample_fields)
        emitter.emit(b)

        for portion in sample.get("portions", []):
            for analyte in portion.get("analytes", []):
                for aliquot in analyte.get("aliquots", []):
                    aliquot_fields = extract(aliquot, ["analyte_type", "submitter_id"])
                    fields = dict(sample_fields)
                    fields.update(aliquot_fields)
                    a = Biosample(aliquot["aliquot_id"], fields)
                    emitter.emit(a)

emitter.close()
