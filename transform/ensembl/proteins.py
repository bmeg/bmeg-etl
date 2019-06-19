import csv

import bmeg.ioutils
from bmeg import Protein, Transcript, Project, Protein_Transcript_Transcript
from bmeg.emitter import JSONEmitter


PROJECT_ID = Project.make_gid("Reference")
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
            p = Protein(submitter_id=Protein.make_gid(protein_id),
                        protein_id=protein_id,
                        uniprot_id=uniprot_id,
                        genome=GENOME_BUILD,
                        project_id=PROJECT_ID)
            emitter.emit_vertex(p)
            emitter.emit_edge(
                Protein_Transcript_Transcript(
                    from_gid=p.gid(),
                    to_gid=Transcript.make_gid(transcript_id)
                ),
                emit_backref=True
            )
            emitted_proteins.append(protein_id)

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
