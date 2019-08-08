"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""

from gdcutils import query_gdc, get_file
from bmeg.utils import ensure_directory
import json
import os


# The GDC API requires you to request that nested fields be expanded.
# https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/#cases-field-groups
#
# Note that (as of this writing) we are expanding most but
# not all possible fields. Mostly we're skipping "files" data.
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


def compounds(output_dir='source/gdc/compounds', parameters={}):
    """ the only way to get drugs is to download files and parse them"""
    my_filters = json.loads("""
    {"op":"and","content":[{"op":"in","content":{"field":"files.data_type","value":["Clinical data"]}},{"op":"in","content":{"field":"files.tags","value":["drug"]}}]}
    """)
    if 'filters' in parameters:
        original_content = parameters['filters']
        my_filters['content'].append(original_content)
    parameters['filters'] = my_filters

    ensure_directory(output_dir)

    for row in query_gdc("legacy/files", parameters):
        get_file(row['file_id'], '{}/{}.tsv'.format(output_dir, row['file_id']))


def cases(output_file='source/gdc/cases.json', parameters={}):
    # Crawl all cases, samples, aliquots to generate
    # BMEG Cases, Samples, and Aliquots.
    ensure_directory(os.path.dirname(output_file))
    fh = open(output_file, 'w')
    for row in query_gdc("cases", {"expand": expand_case_fields}):
        fh.write(json.dumps(row))
        fh.write("\n")


if __name__ == "__main__":
    cases()
    compounds()
