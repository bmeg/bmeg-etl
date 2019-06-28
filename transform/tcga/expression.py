from collections import defaultdict
import csv
from glob import glob
import gzip
import os.path
import logging
import subprocess
import sys

from bmeg import (GeneExpression, TranscriptExpression, Aliquot, Project,
                  GeneExpression_Aliquot_Aliquot, TranscriptExpression_Aliquot_Aliquot)

from bmeg.emitter import JSONEmitter
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging


def transform(source_path,
              id_map_file="source/tcga/expression/transcript-level/TCGA_ID_MAP.csv",
              gene_map_file="source/ensembl/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.trans_gene.tsv",
              emitter_directory="tcga"):

    # check if we are doing one file at time
    p, file_name = os.path.split(source_path)
    prefix = file_name.split('_')[1]
    emitter = JSONEmitter(directory=emitter_directory, prefix=prefix)
    logging.debug('individual file prefix {}'.format(prefix))

    # Map CGHub analysis IDs to GDC Aliquot IDs
    r = csv.DictReader(open(id_map_file))
    id_map = {}
    project_map = {}
    for row in r:
        id_map[row["CGHubAnalysisID"]] = row["Aliquot_id"]
        project_map[row["CGHubAnalysisID"]] = "TCGA-" + row["Disease"]

    r = csv.reader(open(gene_map_file), delimiter="\t")
    gene_map = {row[0]: row[1] for row in r}

    reader = csv.reader(gzip.open(source_path, "rt"), delimiter="\t")
    header = next(reader)
    samples = header[1:]

    # collect expression for all aliquots and transcripts
    collect = defaultdict(dict)

    for row in reader:
        feature_ids = row[0].split("|")
        transcript_id = feature_ids[0].split(".")[0]

        for cghub_id, raw_expr in zip(samples, row[1:]):
            expr = float(raw_expr)
            aliquot_id = id_map[cghub_id]
            project_id = project_map[cghub_id]
            collect[aliquot_id][transcript_id] = expr

    for aliquot_id, values in collect.items():
        t = TranscriptExpression(
            id=TranscriptExpression.make_gid(aliquot_id),
            metric="TPM",
            method="Illumina Hiseq",
            values=values,
            project_id=Project.make_gid(project_id)
        )
        emitter.emit_vertex(t)
        emitter.emit_edge(
            TranscriptExpression_Aliquot_Aliquot(
                from_gid=t.gid(),
                to_gid=Aliquot.make_gid(aliquot_id)
            ),
            emit_backref=True
        )

        geneValues = {}
        for k, v in values.items():
            if k not in gene_map:
                logging.info("%s=%f not found in mapping" % (k, v))
            else:
                gene = gene_map[k]
                geneValues[gene] = geneValues.get(gene, 0) + v

        g = GeneExpression(
            id=GeneExpression.make_gid(aliquot_id),
            metric="GENE_TPM",
            method="Illumina Hiseq",
            values=geneValues,
            project_id=Project.make_gid(project_id)
        )
        emitter.emit_vertex(g)
        emitter.emit_edge(
            GeneExpression_Aliquot_Aliquot(
                from_gid=g.gid(),
                to_gid=Aliquot.make_gid(aliquot_id)
            ),
            emit_backref=True
        )

    emitter.close()


def make_parallel_workstream(source_path, jobs, dry_run=False):
    """ equivalent of: ls -1 source/tcga/expression/transcript-level/TCGA_*_tpm.tsv.gz | parallel --jobs 10 python3.7 transform/tcga/expression.py --source_path """
    with open('/tmp/tcga_expression_transform.txt', 'w') as outfile:
        for path in glob(source_path):
            outfile.write("{}\n".format(path))
    try:
        cmd = 'cat /tmp/tcga_expression_transform.txt | parallel --jobs {} {} {} --source_path'.format(jobs, sys.executable, __file__)
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
    parser.add_argument("--source_path", default="source/tcga/expression/transcript-level/*_tpm.tsv.gz", help="path to file(s)")
    parser.add_argument("--jobs", default=2, help="number of jobs to run in parallel")
    options = parser.parse_args()
    default_logging(options.loglevel)
    if '*' in options.source_path:
        make_parallel_workstream(source_path=options.source_path, jobs=options.jobs)
    else:
        transform(source_path=options.source_path)
