import ijson

from bmeg import (Protein, PfamFamily, ProteinStructure, Project,
                  Protein_PfamFamilies_PfamFamily,
                  Protein_ProteinStructures_ProteinStructure)

from bmeg.emitter import JSONEmitter


# AS-IS: ftp://ftp.ensembl.org/pub/release-93/json/homo_sapiens/homo_sapiens.json
# TODO ?: ftp://ftp.ensembl.org/pub/grch37/release-94/rdf/homo_sapiens.ttl.gz
def transform(data_path='source/pfam/homo_sapiens.json',
              emitter_prefix=None,
              emitter_directory='pfam'):

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
                            for feature in translation.get('protein_features', []):
                                if feature.get('dbname', "") == "Pfam":
                                    if (prot_id, feature['name']) not in dedup:
                                        emitter.emit_edge(
                                            Protein_PfamFamilies_PfamFamily(
                                                from_gid=prot_id,
                                                to_gid=PfamFamily.make_gid(feature['name']),
                                                data={'start': int(feature['start']), 'end': int(feature['end'])}
                                            ),
                                            emit_backref=True
                                        )
                                        dedup[(prot_id, feature['name'])] = None

                            for pdb in translation.get('PDB', []):
                                if (prot_id, pdb) not in dedup:
                                    emitter.emit_edge(
                                        Protein_ProteinStructures_ProteinStructure(
                                            from_gid=prot_id,
                                            to_gid=ProteinStructure.make_gid(pdb),
                                        ),
                                        emit_backref=True
                                    )
                                    dedup[(prot_id, pdb)] = None
                                if pdb not in dedup:
                                    emitter.emit_vertex(
                                        ProteinStructure(
                                            id=ProteinStructure.make_gid(pdb),
                                            structure_id=pdb,
                                            project_id=Project.make_gid("Reference")
                                        )
                                    )
                                    dedup[pdb] = None
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    transform()
