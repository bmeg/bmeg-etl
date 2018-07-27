"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""

import json

import bmeg.models.emitter
from bmeg.models.vertex_models import Individual, Biosample
from bmeg.requests import Client

URL_BASE = "https://api.gdc.cancer.gov/"
client = Client()
#emitter = JSONEmitter("gdc")
#emitter = BSONEmitter("gdc")
#emitter = MsgpackEmitter("gdc")
emitter = bmeg.models.emitter.DebugEmitter()

# The GDC API requires you to explicitly request that nested fields be expanded.
# https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/#cases-field-groups
#
# Note that (as of this writing) we are expanding most but not all possible fields.
# Mostly we're skipping "files" data.
expand = """
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
""".strip().split()

# These are the fields we want to keep from the GDC Case (BMEG Individual).
GDC_CASE_FIELDS = """
diagnoses
demographic
disease_type
primary_site
summary
""".strip().split()

def query_gdc(query, params):
    """
    query_gdc makes a query to the GDC API while handling common issues
    like pagination, retries, etc.

    The return value is an iterator.
    """
    # Copy input params to avoid modification.
    params = dict(params)
    page_size = 100
    params['size'] = page_size

    # Iterate through all the pages.
    while True:
        req = client.get(query, params=params)
        data = req.json()['data']
        # print(json.dumps(data, indent=4))

        hits = data.get("hits", [])
        if len(hits) == 0:
            return

        for hit in hits:
            yield hit

        # Get the next page.
        params['from'] = data['pagination']['from'] + page_size


def extract(data, keys):
    return {key: data.get(key) for key in keys}

for row in query_gdc(URL_BASE + "cases", {"expand": ",".join(expand)}):
    #print(json.dumps(row, indent=True))

    i = Individual(row["id"], extract(row, GDC_CASE_FIELDS))
    emitter.emit(i)

    for sample in row.get("samples", []):

        sample_fields = extract(row, ["tumor_descriptor", "sample_type", "submitter_id"])
        b = Biosample(sample["sample_id"], sample_fields)
        emitter.emit(b)

        for portion in sample.get("portions", []):
            for analyte in portion.get("analytes", []):
                for aliquot in portion.get("aliquots", []):
                    aliquot_fields = extract(row, ["analyte_type", "submitter_id"])
                    fields = dict(sample_fields)
                    fields.update(aliquot_fields)
                    a = Biosample(aliquot["aliquot_id"], fields)
                    emitter.emit(a)

emitter.close()
