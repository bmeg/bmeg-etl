import csv
import gzip
import json

from bmeg.vertex import RPPA, Aliquot
from bmeg.edge import RPPAOf
from bmeg.emitter import JSONEmitter
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging


def transform(source_path,
              aliquot_path="outputs/gdc/Aliquot.Vertex.json.gz",
              emitter_directory="rppa"):

    aliquot_lookup = {}
    with gzip.open(aliquot_path) as handle:
        for line in handle:
            a = json.loads(line)
            barcode = a["data"]["gdc_attributes"]["submitter_id"]
            gdc_id = a["data"]["gdc_attributes"]["aliquot_id"]
            aliquot_lookup[barcode] = gdc_id

    emitter = JSONEmitter(directory=emitter_directory, prefix=None)
    reader = csv.reader(open(source_path, "rt"), delimiter=",")
    header = next(reader)
    feature_ids = header[2:]
    mapped_features = []
    for f in feature_ids:
        # TODO need to map antibody ids to proteins
        # Below goes from antibody --> gene symbol...
        # https://doi.org/10.1371/journal.pone.0188016.s004 --> journal.pone.0188016.s004.xlsx
        gene_symbol = antibody2genesymbol[f]
        mapped_features.append(protein_id)

    for row in reader:
        aliquot_id = aliquot_lookup[row[0]]
        vals = zip(mapped_features, row[2:])
        r = RPPA(
            id=aliquot_id,
            source="tcpa",
            method="",
            values=vals
        )
        emitter.emit_vertex(r)
        emitter.emit_edge(
            RPPAOf(),
            from_gid=r.gid(),
            to_gid=Aliquot.make_gid(aliquot_id)
        )

    emitter.close()


if __name__ == "__main__":
    parser = default_argument_parser()
    parser.add_argument("--source_path", default="source/rppa/TCGA-PANCAN32-L4.csv", help="path to source file")
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform(source_path=options.source_path)
