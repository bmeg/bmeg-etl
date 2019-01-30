import ijson

from bmeg.vertex import Protein, PFAMFamily, ProteinStructure
from bmeg.edge import PFAMAlignment, StructureFor
from bmeg.emitter import JSONEmitter


# AS-IS: ftp://ftp.ensembl.org/pub/release-93/json/homo_sapiens/homo_sapiens.json
# TODO ?: ftp://ftp.ensembl.org/pub/grch37/release-94/rdf/homo_sapiens.ttl.gz
def transform(data_path='source/pfam/homo_sapiens.json',
              emitter_directory='pfam'):

    emitter = JSONEmitter(directory=emitter_directory, prefix=None)

    dedup = []
    with open(data_path) as handle:
        genes = ijson.items(handle, 'genes.item')
        for gene in genes:
            for transcript in gene['transcripts']:
                if 'translations' in transcript:
                    for translation in transcript['translations']:
                        if translation['id'].startswith("ENSP"):
                            prot = Protein(protein_id=translation['id'].split(".")[0], transcript_id=transcript['id'].split(".")[0], genome="GRCh38")
                            for feature in translation['protein_features']:
                                if feature.get('dbname', "") == "Pfam":
                                    emitter.emit_edge(
                                        PFAMAlignment(start=int(feature['start']), end=int(feature['end'])),
                                        from_gid=prot.gid(),
                                        to_gid=PFAMFamily.make_gid(feature['name'])
                                    )
                            for pdb in translation.get('PDB', []):
                                emitter.emit_edge(
                                    StructureFor(),
                                    from_gid=ProteinStructure.make_gid(pdb),
                                    to_gid=prot.gid()
                                )
                                if pdb in dedup:
                                    continue
                                emitter.emit_vertex(ProteinStructure(pdb))
                                dedup.append(pdb)
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
