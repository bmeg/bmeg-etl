from collections import defaultdict
import csv
from glob import iglob
import gzip
from types import SimpleNamespace as SN
import ujson

from bmeg.vertex import Expression, Aliquot, ExpressionMetric, Biosample, Individual, Project
from bmeg.edge import ExpressionOf, AliquotFor, InProject, BiosampleFor
from bmeg.emitter import JSONEmitter
import bmeg.ioutils
from bmeg.util.logging import default_logging
from bmeg.util.cli import default_argument_parser


def transform(source_path,
              id_map_file="source/tcga/expression/transcript-level/TCGA_ID_MAP.csv",
              emitter_directory="ccle",
              prefix='tatlow',
              biosample_path='outputs/ccle/Biosample.Vertex.json*'):

    emitter = JSONEmitter(directory=emitter_directory, prefix=prefix)

    reader = csv.reader(gzip.open(source_path, "rt"), delimiter="\t")
    header = next(reader)
    samples = header[1:]

    # lookup table of Broad_ID
    bio_samples = {}
    for path in iglob(biosample_path):  # match .json or .json.gz
        input_stream = bmeg.ioutils.reader(path)
        for line in input_stream:
            biosample = SN(**ujson.loads(line))
            ccle_attributes = SN(**biosample.data['ccle_attributes'])
            # mangle the lookup keys
            bio_samples[ccle_attributes.CCLE_Name.lower()] = ccle_attributes.Broad_ID
            bio_samples[ccle_attributes.CCLE_Name.split('_')[0].lower()] = ccle_attributes.Broad_ID
            bio_samples[ccle_attributes.Aliases.lower()] = ccle_attributes.Broad_ID
        break
    assert len(bio_samples.keys()), 'No biosamples found, does {} exist?'.format(biosample_path)

    # collect expression for all aliquots and transcripts
    collect = defaultdict(dict)
    missing_cell_lines = {}
    for row in reader:
        feature_ids = row[0].split("|")
        transcript_id = feature_ids[0]

        for cghub_id, raw_expr in zip(samples, row[1:]):
            expr = float(raw_expr)
            # match to existing biosample/aliquot
            aliquot_id = None
            # mangle the lookup keys
            ccle_id = cghub_id.split('.')[1]
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
            keys = [k.lower() for k in keys]
            for k in keys:
                if k in bio_samples:
                    aliquot_id = bio_samples[k]
                    break

            # if no match, we will need to create project->individual->biosample->aliquot
            if not aliquot_id:
                aliquot_id = missing_cell_lines.get(ccle_id, None)
                if not aliquot_id:
                    aliquot_id = ccle_id
                    missing_cell_lines[ccle_id] = aliquot_id
            # if we matched, collect it for output
            collect[aliquot_id][transcript_id] = expr

    # render the expression values
    for aliquot_id, values in collect.items():
        g = Expression(
            id=aliquot_id,
            source="ccle",
            metric=ExpressionMetric.TPM,
            method="Illumina Hiseq",
            values=values,
        )
        emitter.emit_vertex(g)
        emitter.emit_edge(
            ExpressionOf(),
            from_gid=g.gid(),
            to_gid=Aliquot.make_gid(aliquot_id)
        )

    # render missing cell lines
    individual_gids = project_gids = []
    for ccle_id, aliquot_id in missing_cell_lines.items():
        b = Biosample(aliquot_id)
        emitter.emit_vertex(b)

        a = Aliquot(aliquot_id=aliquot_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid(),
        )

        i = Individual(individual_id='CCLE:{}'.format(aliquot_id))
        if i.gid() not in individual_gids:
            emitter.emit_vertex(i)
            individual_gids.append(i.gid())
        emitter.emit_edge(
            BiosampleFor(),
            b.gid(),
            i.gid(),
        )

        # first see if we have wholesale name changes
        project_id = ccle_id
        # strip off prefix
        name_parts = project_id.split('_')
        name_start = 1
        if len(name_parts) == 1:
            name_start = 0
        project_id = '_'.join(name_parts[name_start:])
        # create project
        p = Project(project_id='CCLE:{}'.format(project_id))
        if p.gid() not in project_gids:
            emitter.emit_vertex(p)
            project_gids.append(p.gid())
        emitter.emit_edge(
            InProject(),
            i.gid(),
            p.gid(),
        )

    emitter.close()


if __name__ == "__main__":
    parser = default_argument_parser()
    parser.add_argument("--source_path", default="source/ccle/expression/CCLE_tpm.tsv.gz", help="path to file")
    options = parser.parse_args()
    default_logging(options.loglevel)
    transform(source_path=options.source_path)
