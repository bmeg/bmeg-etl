"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""

from bmeg.util.cli import default_argument_parser
from bmeg.edge import InProject, BiosampleFor, AliquotFor
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Individual, Biosample, Project, Aliquot

from transform.gdc.gdcutils import extract, query_gdc


parser = default_argument_parser()

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

# These are the fields we want to keep from the GDC Case (BMEG Individual).
keep_case_fields = """
diagnoses
demographic
disease_type
primary_site
summary
project
""".strip().split()


def transform(emitter, parameters={}):
    # Crawl all cases, samples, aliquots to generate
    # BMEG Individuals, Biosamples, and Aliquots.
    parameters['expand'] = expand_case_fields
    for row in query_gdc("cases", parameters):
        i = Individual(row["id"], extract(row, keep_case_fields))
        emitter.emit_vertex(i)
        emitter.emit_edge(
            InProject(),
            i.gid(),
            Project.make_gid(i.gdc_attributes["project"]["project_id"]),
        )

        for sample in row.get("samples", []):
            sample_fields = extract(
                sample,
                ["tumor_descriptor", "sample_type", "submitter_id"],
            )
            b = Biosample(sample["sample_id"], sample_fields)
            emitter.emit_vertex(b)

            emitter.emit_edge(
                BiosampleFor(),
                b.gid(),
                i.gid(),
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
                        a = Aliquot(aliquot_id=aliquot["aliquot_id"], gdc_attributes=fields)
                        emitter.emit_vertex(a)

                        emitter.emit_edge(
                            AliquotFor(),
                            a.gid(),
                            b.gid(),
                        )


if __name__ == "__main__":
    args = parser.parse_args()
    emitter = JSONEmitter("gdc")
    transform(emitter)
    emitter.close()
