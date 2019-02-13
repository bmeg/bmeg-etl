"""
Bulk download case, sample, project, etc. data from GDC.
https://gdc.cancer.gov/
"""

from bmeg.util.cli import default_argument_parser
from bmeg.edge import InProject, SampleFor, AliquotFor, TreatedWith
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Case, Sample, Project, Aliquot
from bmeg.ioutils import read_tsv
from transform.gdc.gdcutils import extract, query_gdc, get_file
from bmeg.enrichers.drug_enricher import compound_factory
import json
import logging
import os


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

# These are the fields we want to keep from the GDC Case (BMEG Case).
keep_case_fields = """
diagnoses
demographic
disease_type
primary_site
summary
project
""".strip().split()


def compounds(emitter, parameters={}, output_path='/tmp', case_gids=None):
    """ the only way to get drugs is to download files and parse them"""
    my_filters = json.loads("""
    {"op":"and","content":[{"op":"in","content":{"field":"files.data_type","value":["Clinical data"]}},{"op":"in","content":{"field":"files.tags","value":["drug"]}}]}
    """)
    if 'filters' in parameters:
        original_content = parameters['filters']
        my_filters['content'].append(original_content)
    parameters['filters'] = my_filters
    logging.warning(parameters)
    for row in query_gdc("legacy/files", parameters):
        file_id = row['file_id']
        path = get_file(file_id, '{}/{}'.format(output_path, file_id))
        tsv_in = read_tsv(path)
        c = 0
        dedupe = []
        dedupe_compounds = []
        for line in tsv_in:
            c += 1
            if c < 3:  # skip three header lines
                continue
            bcr_patient_uuid = line['bcr_patient_uuid']
            pharmaceutical_therapy_drug_name = line['pharmaceutical_therapy_drug_name']
            case_gid = Case.make_gid(bcr_patient_uuid)
            # is this a requested case?
            if case_gids is None or case_gid in case_gids:
                t = (bcr_patient_uuid, pharmaceutical_therapy_drug_name)
                # have we seen this combination before?
                if t not in dedupe:
                    compound = compound_factory(name=pharmaceutical_therapy_drug_name)
                    compound_gid = compound.gid()
                    # have we seen this compound before?
                    if compound_gid not in dedupe_compounds:
                        emitter.emit_vertex(compound)
                        dedupe_compounds.append(compound_gid)
                    emitter.emit_edge(
                        TreatedWith(),
                        case_gid,
                        compound_gid,
                    )
                    dedupe.append(t)
        os.remove(path)


def transform(emitter, parameters={}):
    # Crawl all cases, samples, aliquots to generate
    # BMEG Cases, Samples, and Aliquots.
    parameters['expand'] = expand_case_fields
    case_gids = []
    for row in query_gdc("cases", parameters):
        c = Case(row["id"], extract(row, keep_case_fields))
        emitter.emit_vertex(c)
        case_gid = c.gid()
        case_gids.append(case_gid)
        emitter.emit_edge(
            InProject(),
            case_gid,
            Project.make_gid(c.gdc_attributes["project"]["project_id"]),
        )

        for sample in row.get("samples", []):
            sample_fields = extract(
                sample,
                ["tumor_descriptor", "sample_type", "submitter_id"],
            )
            s = Sample(sample["sample_id"], sample_fields)
            emitter.emit_vertex(s)

            emitter.emit_edge(
                SampleFor(),
                s.gid(),
                c.gid(),
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
                            s.gid(),
                        )
    # now use the file endpoint to get compounds
    compounds(emitter, parameters, case_gids=case_gids)


if __name__ == "__main__":
    args = parser.parse_args()
    emitter = JSONEmitter(directory="gdc")
    transform(emitter)
    emitter.close()
