import ijson

from bmeg import (Protein, PfamFamily, ProteinStructure,
                  Protein_PfamFamilies_PfamFamily,
                  Protein_ProteinStructures_ProteinStructure)

from bmeg.emitter import JSONEmitter


# AS-IS: ftp://ftp.ensembl.org/pub/release-93/json/homo_sapiens/homo_sapiens.json
# TODO ?: ftp://ftp.ensembl.org/pub/grch37/release-94/rdf/homo_sapiens.ttl.gz
def transform(data_path='source/pfam/homo_sapiens.json',
              emitter_directory='pfam',
              emitter_prefix='toproteins'):

    emitter = JSONEmitter(directory=emitter_directory, prefix=emitter_prefix)

    dedup = {}
    with open(data_path) as handle:
        genes = ijson.items(handle, 'genes.item')
        for gene in genes:
            for transcript in gene['transcripts']:
                if 'translations' in transcript:
                    for translation in transcript['translations']:
                        if translation['id'].startswith("ENSP"):
                            prot_id = Protein.make_gid(translation['id'].split(".")[0])
                            for feature in translation['protein_features']:
                                if feature.get('dbname', "") == "Pfam":
                                    emitter.emit_edge(
                                        Protein_PfamFamilies_PfamFamily(
                                            from_gid=prot_id,
                                            to_gid=PfamFamily.make_gid(feature['name']),
                                            data={'start': int(feature['start']), 'end': int(feature['end'])}
                                        ),
                                        emit_backref=True
                                    )

                            for pdb in translation.get('PDB', []):
                                emitter.emit_edge(
                                    Protein_ProteinStructures_ProteinStructure(
                                        from_gid=prot_id,
                                        to_gid=ProteinStructure.make_gid(pdb),
                                    ),
                                    emit_backref=True
                                )
                                if pdb in dedup:
                                    continue
                                emitter.emit_vertex(
                                    ProteinStructure(
                                        submitter_id=ProteinStructure.make_gid(pdb),
                                        structure_id=pdb,
                                        project_id="Reference"
                                    )
                                )
                                dedup[pdb] = None
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
