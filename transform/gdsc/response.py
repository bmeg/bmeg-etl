import re

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Aliquot, DrugResponse, DrugResponseMetric
from bmeg.edge import DrugResponseIn, ResponseTo
from bmeg.util.logging import default_logging
from bmeg.util.cli import default_argument_parser
from bmeg.enrichers.drug_enricher import compound_factory
import logging

DEFAULT_PREFIX = 'gdsc'
DEFAULT_DIRECTORY = 'gdsc'

BROAD_LOOKUP = {
    # missing broad ids
    "A253_UPPER_AERODIGESTIVE_TRACT": "ACH-000740",
    "RPMI2650_UPPER_AERODIGESTIVE_TRACT": "ACH-001385",
    "U698M_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE": "ACH-001680",
    # incorrect broad ids
    # GDSC
    # ACH-001311,HSC2_UPPER_AERODIGESTIVE_TRACT
    # ACH-001312,HSC3_UPPER_AERODIGESTIVE_TRACT
    # ACH-001319,SCC25_UPPER_AERODIGESTIVE_TRACT
    # ACH-001320,SCC4_UPPER_AERODIGESTIVE_TRACT
    # ACH-001821,EFM19_BREAST

    "HSC2_UPPER_AERODIGESTIVE_TRACT": "ACH-000472",
    "HSC3_UPPER_AERODIGESTIVE_TRACT": "ACH-000778",
    "SCC25_UPPER_AERODIGESTIVE_TRACT": "ACH-000188",
    "SCC4_UPPER_AERODIGESTIVE_TRACT": "ACH-000238",
    "EFM19_BREAST": "ACH-000330",
}


def transform(
        path="source/gdsc/GDSC_AUC.csv",
        emitter_prefix=DEFAULT_PREFIX,
        emitter_directory=DEFAULT_DIRECTORY,
):
    logging.info('transform')
    emitter = JSONEmitter(prefix=emitter_prefix, directory=emitter_directory)
    r = bmeg.ioutils.read_csv(path)

    # Fix up the headers:

    # Depmap added some strange headers and changed their sample ID codes,
    # The GDSC data currently is downloaded from Depmap.org.
    # Depmap started a new sample ID type (Broad ID) in order to ensure uniqueness.
    # We're not using Broad IDs yet, so we parse the CCLE sample ID
    # out of the header, and reset the csv reader's fieldnames.
    rx = re.compile("^(.*) \((.*)\)$")

    # The first column header is blank.
    replace_with = ["compound_name"]

    for field in r.fieldnames[1:]:
        m = rx.search(field)
        sample_code = None
        if not m:
            # Some don't have a Broad ID.
            sample_code = BROAD_LOOKUP.get(field, None)
        else:
            sample_code = BROAD_LOOKUP.get(m.group(1), m.group(2))

        assert sample_code, 'no broad id for {}'.format(field)
        replace_with.append(sample_code)

    assert len(replace_with) == len(r.fieldnames)
    r.fieldnames = replace_with

    # Iterate all rows, writing out the expression for different tissue types
    # to separate files.
    c = t = 0
    compound_gids = []
    for row in r:
        c += 1
        t += 1
        for key, raw_value in row.items():
            # Skip the first column
            if key == "compound_name":
                continue

            sample = key
            value = None
            if raw_value != "NA":
                value = float(raw_value)

            e = DrugResponse(
                compound_name=row["compound_name"],
                sample_id=sample,
                metric=DrugResponseMetric.AUC,
                value=value,
            )
            emitter.emit_vertex(e)
            emitter.emit_edge(
                DrugResponseIn(),
                e.gid(),
                Aliquot.make_gid(sample),
            )
            # create compound
            compound = compound_factory(name=row["compound_name"])
            if compound.gid() not in compound_gids:
                emitter.emit_vertex(compound)
                compound_gids.append(compound.gid())

            # and an edge to it
            emitter.emit_edge(
                ResponseTo(),
                e.gid(),
                compound.gid(),
            )

            if c % 10 == 0:
                logging.info('imported {}'.format(t))
                c = 0

    emitter.close()


if __name__ == "__main__":
    parser = default_argument_parser()
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform()
