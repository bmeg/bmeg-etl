import re

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample, DrugResponse, DrugResponseMetric
from bmeg.edge import DrugResponseIn, ResponseTo
from bmeg.util.logging import default_logging
from bmeg.util.cli import default_argument_parser
from bmeg.enrichers.drug_enricher import enrich
import logging

DEFAULT_PREFIX = 'gdsc'


def transform(path="source/gdsc/GDSC_AUC.csv", prefix=DEFAULT_PREFIX):
    logging.info('transform')
    emitter = JSONEmitter(prefix)
    r = bmeg.ioutils.read_csv(path)

    # Fix up the headers:

    # Depmap added some strange headers and changed their sample ID codes,
    # The GDSC data currently is downloaded from Depmap.org.
    # Depmap started a new sample ID type (Broad ID) in order to ensure uniqueness.
    # We're not using Broad IDs yet, so we parse the CCLE sample ID
    # out of the header, and reset the csv reader's fieldnames.
    rx = re.compile("^(.*) \((.*)\)$")  # noqa W605 invalid escape sequence

    # The first column header is blank.
    replace_with = ["compound_name"]

    for field in r.fieldnames[1:]:
        m = rx.search(field)
        # Some don't have a Broad ID.
        if not m:
            replace_with.append(field)
            continue

        sample_code = m.group(1)
        replace_with.append(sample_code)

    assert len(replace_with) == len(r.fieldnames)
    r.fieldnames = replace_with

    # Iterate all rows, writing out the expression for different tissue types
    # to separate files.
    c = t = 0
    compound_cache = {}
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
                Biosample.make_gid(sample),
            )
            # create compound
            compound = compound_cache.get(row["compound_name"], None)
            if not compound:
                compound = enrich(row["compound_name"])
                emitter.emit_vertex(compound)
                compound_cache[row["compound_name"]] = compound
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
