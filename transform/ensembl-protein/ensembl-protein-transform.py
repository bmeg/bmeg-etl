#!/usr/bin/env python

import ijson

from bmeg.vertex import Protein, PFAMFamily, Transcript, ProteinStructure
from bmeg.edge import PFAMAlignment, ProteinFor, StructureFor
from bmeg.emitter import JSONEmitter

# AS-IS: ftp://ftp.ensembl.org/pub/release-93/json/homo_sapiens/homo_sapiens.json
# TODO ?: ftp://ftp.ensembl.org/pub/grch37/release-94/rdf/homo_sapiens.ttl.gz

emitter = JSONEmitter("ensembl-protein")

path = "source/ensembl-protein/homo_sapiens.json"
dedup = []
with open(path) as handle:
    genes = ijson.items(handle, 'genes.item')
    for gene in genes:
        for transcript in gene['transcripts']:
            if 'translations' in transcript:
                for translation in transcript['translations']:
                    if translation['id'].startswith("ENSP"):
                        prot = Protein(protein_id=translation['id'], transcript_id=transcript['id'])
                        emitter.emit_vertex(prot)
                        emitter.emit_edge(
                            ProteinFor(),
                            from_gid=prot.gid(),
                            to_gid=Transcript.make_gid(transcript['id'])
                        )
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
