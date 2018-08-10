import bmeg.ioutils
from bmeg.util.cli import default_argument_parser
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample, Project


parser = default_argument_parser()
parser.add_argument("samples_file")


def transform(emitter, reader):
    emitter.emit_vertex(Project("gtex", gdc_attributes={}))

    for row in reader:
        b = Biosample(
            row["SAMPID"],
            gdc_attributes={},
            gtex_attributes=row,
        )
        emitter.emit_vertex(b)


if __name__ == "__main__":
    args = parser.parse_args()
    emitter = JSONEmitter("gtex")
    transform(emitter, bmeg.ioutils.read_tsv(args.samples_file))
    emitter.close()
