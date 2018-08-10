import csv
import re

import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.util import cli
from bmeg.vertex import Biosample, DrugResponse, DrugResponseMetric
from bmeg.edge import DrugResponseIn


def transform(emitter, auc_file, sample_info_file):
    r = csv.DictReader(auc_file)
    # The first column header is blank
    r.fieldnames[0] = "compound_name"

    # Depmap added some strange headers and changed their sample ID codes,
    # so extract the part we need from the column headers and build a mapping
    # from sample ID (Broad ID) to column header.
    column_to_sample = {}
    rx = re.compile("^.* \((.*)\)$")
    for field in r.fieldnames:
        m = rx.search(field)
        if m:
            sample_code = m.group(1)
            column_to_sample[field] = sample_code

    # Load sample metadata.
    for row in csv.DictReader(sample_info_file):
        sample_id = row["Broad_ID"]
        b = Biosample(sample_id, ccle_attributes=row)
        emitter.emit_vertex(b)

    # Iterate all rows, writing out the expression for different tissue types
    # to separate files.
    for row in r:
        for key, raw_value in row.items():
            if key == "compound_name":
                continue

            try:
                sample = column_to_sample[key]
            except KeyError:
                continue

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


parser = cli.default_argument_parser()
parser.add_argument("auc_file")
parser.add_argument("sample_info_file")

if __name__ == "__main__":
    args = parser.parse_args()
    emitter = JSONEmitter("ccle")
    transform(
        emitter,
        bmeg.ioutils.reader(args.auc_file),
        bmeg.ioutils.reader(args.sample_info_file),
    )
    emitter.close()
