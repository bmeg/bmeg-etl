import csv

import bmeg.ioutils
from bmeg import Protein, Uniprot, Project, Protein_Uniprot_Uniprot
from bmeg.emitter import JSONEmitter


PROJECT_ID = Project.make_gid("Reference")
GENOME_BUILD = "GRCh37"
DEFAULT_DIRECTORY = "ensembl"


def transform(protein_table_path='source/ensembl/Homo_sapiens.GRCh37.85.uniprot.tsv.gz',
              emitter_directory=DEFAULT_DIRECTORY):

    emitter = JSONEmitter(directory=emitter_directory, prefix=None)
    inhandle = bmeg.ioutils.reader(protein_table_path)
    reader = csv.DictReader(inhandle, delimiter="\t")
    emitted = {}
    for line in reader:
        protein_id = line['protein_stable_id']
        uniprot_id = line['xref']

        if uniprot_id == "-" or protein_id == "-":
            continue

        if line['info_type'] == 'DIRECT':
            line['source_identity'] = 100
        elif line['info_type'] == 'SEQUENCE_MATCH':
            line['source_identity'] = int(line['source_identity'])
        elif line['info_type'] == 'MISC':
            continue
        else:
            continue

        if line['source_identity'] < 50:
            continue

        if uniprot_id not in emitted:
            p = Uniprot(id=Uniprot.make_gid(uniprot_id),
                        genome=GENOME_BUILD,
                        project_id=PROJECT_ID)
            emitter.emit_vertex(p)
            emitter.emit_edge(
                Protein_Uniprot_Uniprot(
                    from_gid=Protein.make_gid(protein_id),
                    to_gid=p.gid(),
                    data={'source': line['db_name'],
                          'type': line['info_type'],
                          'sequence_match_percent': int(line['source_identity'])}
                ),
                emit_backref=True
            )
            emitted[uniprot_id] = None

    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
