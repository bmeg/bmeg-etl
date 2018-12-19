from collections import defaultdict
import gzip
import logging
import json
import os.path
import tarfile

import bmeg.enrichers.gene_enricher as gene_enricher
from bmeg.vertex import Methylation, MethylationProbe, Aliquot, Gene
from bmeg.edge import MethylationOf, MethylationProbeFor
from bmeg.emitter import JSONEmitter
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging


def transform(source_path,
              aliquot_path="outputs/gdc/Aliquot.Vertex.json.gz",
              kvstore_path="outputs/tcga/methylation.db",
              emitter_directory="tcga"):

    logging.debug("building aliquot lookup...")
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

    probe_id_idx = 0
    bval_idx = 1
    gene_symbol_idx = 2
    chromosome_idx = 3
    coordinate_idx = 4

    logging.debug("opening source tarball...")
    tar = tarfile.open(source_path, "r|gz")
    i = 0
    for member in tar:
        if not member.isfile() or "HumanMethylation450" not in member.name:
            continue
        logging.debug("processing: %s", member.name)
        f = tar.extractfile(member)
        if not f:
            logging.debug("tar.extractfile returned None for: %s", member.name)
            continue
        aliquot_barcode = f.readline().decode('ascii').strip().split("\t")[1]
        f.readline()

        # collect expression for all aliquots and transcripts
        collect = defaultdict(dict)
        for row in f:
            row = row.decode('ascii').strip().split("\t")
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

            # only emit probe vertices when processing the first file
            if i == 0:
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

            # http://gdac.broadinstitute.org/runs/sampleReports/latest/SKCM_Notifications.html
            if aliquot_barcode == "TCGA-XV-AB01-01A-12D-A408-05":
                aliquot_barcode = "TCGA-XV-AB01-06A-12D-A408-05"
            aliquot_id = aliquot_lookup[aliquot_barcode]
            beta_val = row[bval_idx]
            if beta_val == "NA":
                bval = None
            else:
                bval = float(beta_val)
            collect[aliquot_id][probe_id] = bval

        for aliquot_id, values in collect.items():
            m = Methylation(
                id=aliquot_id,
                source="tcga",
                metric="Methylation beta value",
                method=prefix,
                values=values,
            )
            emitter.emit_vertex(m)
            emitter.emit_edge(
                MethylationOf(),
                from_gid=m.gid(),
                to_gid=Aliquot.make_gid(aliquot_id)
            )

        logging.debug("finished processing: %s", member.name)
        i += 1

    tar.close()
    emitter.close()


if __name__ == "__main__":
    parser = default_argument_parser()
    parser.add_argument("--source_path", default="source/tcga/methylation/IlluminaHumanMethylation450.tar.gz", help="path to methylation tarbar")
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform(source_path=options.source_path)
