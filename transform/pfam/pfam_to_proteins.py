import ijson

from bmeg import (Protein, PFAMFamily, ProteinStructure,
                  Protein_PFAMFamilies_PFAMFamily, PFAMFamily_Proteins_Protein,
                  Protein_ProteinStructure_ProteinStructure, ProteinStructure_Protein_Protein)

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
                                    # TODO: edge features
                                    # PFAMAlignment(start=int(feature['start']), end=int(feature['end'])),
                                    emitter.emit_edge(
                                        Protein_PFAMFamilies_PFAMFamily(
                                            from_gid=prot.gid(),
                                            to_gid=PFAMFamily.make_gid(feature['name'])
                                        )
                                    )
                                    emitter.emit_edge(
                                        PFAMFamily_Proteins_Protein(
                                            from_gid=PFAMFamily.make_gid(feature['name']),
                                            to_gid=prot.gid()
                                        )
                                    )
                            for pdb in translation.get('PDB', []):
                                emitter.emit_edge(
                                    ProteinStructure_Protein_Protein(
                                        from_gid=ProteinStructure.make_gid(pdb),
                                        to_gid=prot.gid()
                                    )
                                )
                                emitter.emit_edge(
                                    Protein_ProteinStructure_ProteinStructure(
                                        from_gid=prot.gid(),
                                        to_gid=ProteinStructure.make_gid(pdb),
                                    )
                                )
                                if pdb in dedup:
                                    continue
                                emitter.emit_vertex(ProteinStructure(pdb))
                                dedup.append(pdb)
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
