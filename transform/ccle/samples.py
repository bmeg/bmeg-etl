import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample, Aliquot
from bmeg.edge import AliquotFor


def transform(path="source/ccle/CCLE_sample_info_file_2012-10-18.txt",
              prefix="ccle"):
    emitter = JSONEmitter(prefix)
    reader = bmeg.ioutils.read_tsv(path)

    # Load sample metadata.
    for row in reader:
        sample_id = row["CCLE name"]
        b = Biosample(sample_id, ccle_attributes=row)
        emitter.emit_vertex(b)
        a = Aliquot(aliquot_id=sample_id)
        emitter.emit_vertex(a)
        emitter.emit_edge(
            AliquotFor(),
            a.gid(),
            b.gid(),
        )
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
