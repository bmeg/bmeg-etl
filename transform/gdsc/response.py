import csv
import re

import bmeg.ioutils
from bmeg.util import logging
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample, DrugResponse, DrugResponseMetric
from bmeg.edge import DrugResponseIn


logging.simple_logging()

emitter = JSONEmitter("gdsc")
r = bmeg.ioutils.read_csv("source/gdsc/GDSC_AUC.csv")

# Fix up the headers:

# Depmap added some strange headers and changed their sample ID codes,
# The GDSC data currently is downloaded from Depmap.org.
# Depmap started a new sample ID type (Broad ID) in order to ensure uniqueness.
# We're not using Broad IDs yet, so we parse the CCLE sample ID out of the header,
# and reset the csv reader's fieldnames.
rx = re.compile("^(.*) \((.*)\)$")

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
for row in r:
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

emitter.close()
