import csv
import gzip
import logging
import json
import os.path

import bmeg.enrichers.gene_enricher as gene_enricher
from bmeg.vertex import Methylation, MethylationProbe, Aliquot, Gene
from bmeg.edge import MethylationOf, MethylationProbeFor
from bmeg.emitter import JSONEmitter
from bmeg.stores import KeyValueStore
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging


def transform(source_path,
              aliquot_path="outputs/gdc/Aliquot.Vertex.json.gz",
              kvstore_path="outputs/tcga/methylation.db",
              emitter_directory="tcga"):

    aliquot_lookup = {}
    with gzip.open(aliquot_path) as handle:
        for line in handle:
            a = json.loads(line)
            barcode = a["data"]["gdc_attributes"]["submitter_id"]
            gdc_id = a["data"]["gdc_attributes"]["aliquot_id"]
            aliquot_lookup[barcode] = gdc_id

    # check if we are doing one file at time
    p, file_name = os.path.split(source_path)
    prefix = file_name.split(".")[0]
    emitter = JSONEmitter(directory=emitter_directory, prefix=prefix)

    reader = csv.reader(gzip.open(source_path, "rt"), delimiter="\t")
    header = next(reader)
    samples = header[1:-3]

    probe_id_idx = 0
    gene_symbol_idx = -3
    chromosome_idx = -2
    coordinate_idx = -1

    # collect methylation for all aliquots and transcripts
    db = KeyValueStore(path=kvstore_path)
    db.index()

    for row in reader:
        probe_id = row[probe_id_idx]
        symbol = row[gene_symbol_idx].split(";")[0]
        if symbol:
            try:
                gene = gene_enricher.get_gene(symbol)
                ensembl_id = gene["ensembl_gene_id"]
                symbol = ensembl_id
            except Exception as e:
                logging.warning(str(e))

        if symbol is None:
            symbol = ""

        p = MethylationProbe(
            id=probe_id,
            target=symbol,
            chromosome=row[chromosome_idx],
            position=int(row[coordinate_idx]),
        )
        emitter.emit_vertex(p)

        if symbol.startswith("ENSG"):
            emitter.emit_edge(
                MethylationProbeFor(),
                from_gid=p.gid(),
                to_gid=Gene.make_gid(symbol)
            )

        for aliquot_barcode, beta_val in zip(samples, row[1:-3]):
            # http://gdac.broadinstitute.org/runs/sampleReports/latest/SKCM_Notifications.html
            if aliquot_barcode == "TCGA-XV-AB01-01A-12D-A408-05":
                aliquot_barcode = "TCGA-XV-AB01-06A-12D-A408-05"
            aliquot_id = aliquot_lookup[aliquot_barcode]
            if beta_val == "NA":
                bval = None
            else:
                bval = float(beta_val)
            vals = db.get(aliquot_id)
            if not vals:
                vals = {}
            vals[probe_id] = bval
            db.put(aliquot_id, vals)

    for aliquot_id in db.all_ids():
        m = Methylation(
            id=aliquot_id,
            source="tcga",
            metric="Methylation beta value",
            method=prefix,
            values=db.get(aliquot_id),
        )
        emitter.emit_vertex(m)
        emitter.emit_edge(
            MethylationOf(),
            from_gid=m.gid(),
            to_gid=Aliquot.make_gid(aliquot_id)
        )

    emitter.close()


if __name__ == "__main__":
    parser = default_argument_parser()
    parser.add_argument("--source_path", default="source/tcga/methylation/IlluminaHumanMethylation450.tsv.gz", help="path to methylation matrix")
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform(source_path=options.source_path)
