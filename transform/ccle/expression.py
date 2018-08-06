from collections import defaultdict
import csv
import re

from bmeg import vertex, emitter
from bmeg.util import cli, logging

def transform(emitter):
    f = open("data/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct")

    # Skip GCT file version comment
    next(f)
    # Skip GCT file shape line
    next(f)

    r = csv.DictReader(f, delimiter="\t")

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
    for row in csv.DictReader(open("data/sample_info.csv")):
        sample_id = row["Broad_ID"]
        b = vertex.Biosample(sample_id, ccle_attributes=row)
        emitter.emit_vertex(b)

    # Iterate all rows, writing out the expression for different tissue types
    # to separate files.
    for row in r:
        for key, value in row.items():
            if key == "gene":
                continue

            try:
                sample = column_to_sample[key]
            except KeyError:
                continue

            e = vertex.Expression(
                gene_id=row["Name"],
                sample_id=sample,
                # TODO enum
                scale="RPKM",
                value=float(value),
            )
            emitter.emit_vertex(e)


parser = cli.default_argument_parser()
parser.add_argument("--emitter", type=str, default="json")

if __name__ == "__main__":
    args = parser.parse_args()
    e = emitter.get_emitter(args.emitter, prefix="ccle")
    transform(e)
    e.close()
