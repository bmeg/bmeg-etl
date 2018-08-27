from collections import defaultdict
import csv
from glob import glob
import os

from bmeg.vertex import Expression, Aliquot, ExpressionMetric
from bmeg.edge import ExpressionOf
from bmeg.emitter import JSONEmitter


emitter = JSONEmitter("tcga")

# Map CGHub analysis IDs to GDC Aliquot IDs
id_map_file = "source/tcga/expression/transcript-level/TCGA_ID_MAP.csv"
r = csv.DictReader(open(id_map_file))
id_map = {row["CGHubAnalysisID"]: row["Aliquot_id"] for row in r}

for path in glob("source/tcga/expression/transcript-level/*_tpm.tsv.gz"):
    reader = csv.reader(open(path), delimiter="\t")
    header = next(reader)
    samples = header[1:]

    # collect expression for all aliquots and transcripts
    collect = defaultdict(dict)

    for row in reader:
        feature_ids = row[0].split("|")
        transcript_id = feature_ids[0]

        for cghub_id, raw_expr in zip(samples, row[1:]):
            expr = float(raw_expr)
            aliquot_id = id_map[cghub_id]
            collect[aliquot_id][transcript_id] = expr

    for aliquot_id, values in collect.items():
        g = Expression(
            id=aliquot_id,
            source="tcga",
            scale=ExpressionMetric.TPM,
            values=values,
        )
        emitter.emit_vertex(g)
        emitter.emit_edge(
            ExpressionOf(),
            from_gid=g.gid(),
            to_gid=Aliquot.make_gid(aliquot_id)
        )

emitter.close()
