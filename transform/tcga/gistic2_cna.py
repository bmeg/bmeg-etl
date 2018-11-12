from collections import defaultdict
import csv
from glob import glob
import os.path
import gzip
import json
import logging
import subprocess
import sys

import bmeg.enrichers.gene_enricher as gene_enricher
from bmeg.vertex import Aliquot, CopyNumberAlteration
from bmeg.edge import CopyNumberAlterationOf
from bmeg.emitter import JSONEmitter
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging


def transform(source_path,
              aliquot_path="outputs/gdc/Aliquot.Vertex.json.gz",
              emitter_directory="tcga",
              source="tcga",
              method="gistic2"):

    # check if we are doing one file at time
    p, file_name = os.path.split(source_path)
    prefix = file_name.split('_')[0]
    emitter = JSONEmitter(directory=emitter_directory, prefix=prefix)
    logging.debug('individual file prefix {}'.format(prefix))

    aliquot_lookup = {}
    with gzip.open(aliquot_path) as handle:
        for line in handle:
            a = json.loads(line)
            barcode = a["data"]["gdc_attributes"]["submitter_id"]
            gdc_id = a["data"]["gdc_attributes"]["aliquot_id"]
            aliquot_lookup[barcode] = gdc_id

    reader = csv.reader(open(source_path, "rt"), delimiter="\t")
    header = next(reader)
    samples = header[3:]

    # collect expression for all aliquots and transcripts
    collect = defaultdict(dict)

    for row in reader:
        feature_ids = row[0].split("|")
        symbol = feature_ids[0]
        try:
            genes = gene_enricher.get_gene(symbol)
            ensembl_id = genes[0]['ensembl_gene_id']
            symbol = ensembl_id
        except Exception as e:
            logging.warning(str(e))
            continue  # do we want to drop these data?

        for aliquot_id, val in zip(samples, row[3:]):
            collect[aliquot_id][symbol] = val

    for aliquot_barcode, values in collect.items():
        aliquot_id = aliquot_lookup[aliquot_barcode]
        cna = CopyNumberAlteration(
            id=aliquot_id,
            source=source,
            method=method,
            values=values,
        )
        emitter.emit_vertex(cna)
        emitter.emit_edge(
            CopyNumberAlterationOf(),
            from_gid=cna.gid(),
            to_gid=Aliquot.make_gid(aliquot_id)
        )

    emitter.close()


def make_parallel_workstream(source_path, jobs, dry_run=False):
    """ equivalent of: ls -1 source/tcga/expression/gistic2-firehose/TCGA-*_all_thresholded.by_genes.txt | parallel --jobs 10 python3.7 transform/tcga/gistic2_cna.py --source_path """
    with open('/tmp/tcga_gistic2_transform.txt', 'w') as outfile:
        for path in glob(source_path):
            outfile.write("{}\n".format(path))
    try:
        cmd = 'cat /tmp/tcga_gistic2_transform.txt | parallel --jobs {} {} {} --source_path'.format(jobs, sys.executable, __file__)
        logging.info('running {}'.format(cmd))
        if not dry_run:
            subprocess.check_output(cmd, shell=True)
        else:
            return cmd
    except subprocess.CalledProcessError as run_error:
        logging.exception(run_error)
        raise ValueError("error code {} {}".format(run_error.returncode, run_error.output))


if __name__ == "__main__":
    parser = default_argument_parser()
    parser.add_argument("--source_path", default="source/tcga/gistic2-firehose/*_all_thresholded.by_genes.txt", help="path to file(s)")
    parser.add_argument("--jobs", default=10, help="number of jobs to run in parallel")
    options = parser.parse_args()
    default_logging(options.loglevel)
    if '*' in options.source_path:
        make_parallel_workstream(source_path=options.source_path, jobs=options.jobs)
    else:
        transform(source_path=options.source_path)
