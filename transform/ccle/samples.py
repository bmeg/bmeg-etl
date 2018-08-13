import bmeg.ioutils
from bmeg.emitter import JSONEmitter
from bmeg.vertex import Biosample


emitter = JSONEmitter("ccle")
path = "source/ccle/CCLE_sample_info_file_2012-10-18.txt"
reader = bmeg.ioutils.read_tsv(path)

# Load sample metadata.
for row in reader:
    sample_id = row["CCLE name"]
    b = Biosample(sample_id, ccle_attributes=row)
    emitter.emit_vertex(b)

emitter.close()
