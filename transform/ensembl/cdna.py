
from bmeg import CDNA
from Bio import SeqIO
import gzip
from bmeg.emitter import JSONEmitter


emitter = JSONEmitter("ensembl", "cdna")

handle = gzip.open("source/ensembl/Homo_sapiens.GRCh37.cdna.all.fa.gz", "rt")
fasta_sequences = SeqIO.parse(handle,'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    transcript = name.split(".")[0]
    if transcript.startswith("ENST"):
        c = CDNA(transcript_id=transcript, sequence=sequence)
        emitter.emit_vertex(c)
