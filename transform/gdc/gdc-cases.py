"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""

from pprint import pprint

from bmeg.models.emitter import Emitter
from bmeg.models.vertex_models import Individual, Biosample
from bmeg.requests import Client

URL_BASE = "https://api.gdc.cancer.gov/"
client = Client()
emitter = Emitter("gdc")

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

for i, row in enumerate(query_gdc(URL_BASE + "cases", {"expand": ",".join(expand)})):
    #pprint(row)
    m = Individual(row["id"], extract(row, ["diagnoses", "demographic", "disease_type", "primary_site", "summary"]))
    emitter.emit(m)

emitter.close()
