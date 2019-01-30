from collections import defaultdict
import csv
import gzip

from bmeg.vertex import TranscriptExpression, Aliquot, ExpressionMetric
from bmeg.edge import TranscriptExpressionOf
from bmeg.emitter import JSONEmitter
from bmeg.util.logging import default_logging
from bmeg.util.cli import default_argument_parser
from bmeg.ccle import build_ccle2depmap_conversion_table, missing_ccle_cellline_factory


def transform(source_path,
              emitter_directory="ccle",
              prefix='tatlow',
              biosample_path='outputs/ccle/Biosample.Vertex.json.gz'):

    emitter = JSONEmitter(directory=emitter_directory, prefix=prefix)

    reader = csv.reader(gzip.open(source_path, "rt"), delimiter="\t")
    header = next(reader)
    samples = header[1:]

    # lookup table of Broad_ID
    bio_samples = build_ccle2depmap_conversion_table(biosample_path)

    # collect expression for all aliquots and transcripts
    collect = defaultdict(dict)
    missing_cell_lines = []
    for row in reader:
        feature_ids = row[0].split("|")
        transcript_id = feature_ids[0]

        for cghub_id, raw_expr in zip(samples, row[1:]):
            expr = float(raw_expr)
            ccle_id = cghub_id.split('.')[1]
            # match to existing biosample/aliquot
            aliquot_id = None
            # mangle the lookup keys
            keys = [
                ccle_id,  # G20460.COR-L24.2  -> COR-L24
                ccle_id.split('_')[0],  # G20460.COR-L24_FOO.2  -> COR-L24
                ccle_id.replace('-', ''),  # G20460.COR-L24.2  -> CORL24
                ccle_id.replace('-', ' '),  # G20460.COR-L24.2  -> CORL24
                ccle_id.replace('_', ' '),  # G20460.COR_L24.2  -> COR L24
                ccle_id.replace('_', ''),  # G20460.COR_L24.2  -> CORL24
            ]
            # missed spelled
            if ccle_id.startswith("Hs_"):
                keys.append('HS{}T'.format(ccle_id.split('_')[1]))
            if ccle_id.startswith("TE_"):
                keys.append('TE{}T'.format(ccle_id.split('_')[1]))
            if ccle_id.startswith("TO_"):
                keys.append('TO{}T'.format(ccle_id.split('_')[1]))
            keys = [k.upper() for k in keys]
            for k in keys:
                if k in bio_samples:
                    aliquot_id = bio_samples[k]
                    break

            # if no match, we will need to create project->individual->biosample->aliquot
            if not aliquot_id:
                aliquot_id = ccle_id
                if ccle_id not in missing_cell_lines:
                    missing_cell_lines.append(ccle_id)

            collect[aliquot_id][transcript_id] = expr

    # render the expression values
    for aliquot_id, values in collect.items():
        g = TranscriptExpression(
            id=aliquot_id,
            source="ccle-tatlow",
            metric=ExpressionMetric.TPM,
            method="Illumina Hiseq",
            values=values,
        )
        emitter.emit_vertex(g)
        emitter.emit_edge(
            TranscriptExpressionOf(),
            from_gid=g.gid(),
            to_gid=Aliquot.make_gid(aliquot_id)
        )

    # generate project, individual, biosample, aliquot for missing cell lines
    missing_ccle_cellline_factory(emitter=emitter,
                                  source="CCLE",
                                  missing_ids=missing_cell_lines,
                                  project_id="CCLE:UNKNOWN")

    emitter.close()


if __name__ == "__main__":
    parser = default_argument_parser()
    parser.add_argument("--source_path", default="source/ccle/expression/CCLE_tpm.tsv.gz", help="path to file")
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform(source_path=options.source_path)
