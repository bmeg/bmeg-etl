import csv

import bmeg.ioutils
from bmeg.vertex import Protein, Transcript
from bmeg.edge import HasProtein
from bmeg.emitter import JSONEmitter

GENOME_BUILD = "GRCh37"
DEFAULT_DIRECTORY = "ensembl"


def transform(protein_table_path='source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz',
              emitter_directory=DEFAULT_DIRECTORY):

    emitter = JSONEmitter(directory=emitter_directory, prefix=None)
    inhandle = bmeg.ioutils.reader(protein_table_path)
    reader = csv.DictReader(inhandle, delimiter="\t")
    emitted_proteins = []
    for line in reader:
        transcript_id = line['transcript_stable_id']
        protein_id = line['protein_stable_id']
        uniprot_id = None
        if line['info_type'] == 'DIRECT':
            uniprot_id = line['xref']
        elif line['info_type'] == 'SEQUENCE_MATCH':
            if line['source_identity'] == 100:
                uniprot_id = line['xref']

        if protein_id != "-" and protein_id not in emitted_proteins:
            p = Protein(protein_id=protein_id,
                        transcript_id=transcript_id,
                        uniprot_id=uniprot_id,
                        genome=GENOME_BUILD)
            emitter.emit_vertex(p)
            emitter.emit_edge(
                HasProtein(),
                to_gid=p.gid(),
                from_gid=Transcript.make_gid(transcript_id)
            )
            emitted_proteins.append(protein_id)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
