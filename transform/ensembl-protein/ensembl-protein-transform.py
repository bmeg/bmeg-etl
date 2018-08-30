#!/usr/bin/env python

import sys
import json
import ijson

from bmeg.vertex import Protein, PFAMFamily, Transcript
from bmeg.edge import PFAMAlignment, ProteinFor
from bmeg.emitter import JSONEmitter

# ftp://ftp.ensembl.org/pub/release-93/json/homo_sapiens/homo_sapiens.json

emitter = JSONEmitter("ensembl-protein")

with open(sys.argv[1]) as handle:
    genes = ijson.items(handle, 'genes.item')
    for gene in genes:
        #print(json.dumps(gene, indent=4))
        for transcript in gene['transcripts']:
            if 'translations' in transcript:
                for translation in transcript['translations']:
                    if translation['id'].startswith("ENSP"):
                        prot = Protein(protein_id=translation['id'], transcript_id=transcript['id'])
                        emitter.emit_vertex(prot)
                        emitter.emit_edge(ProteinFor(),
                            from_gid=prot.gid(),
                            to_gid=Transcript.make_gid(transcript['id'])
                        )
                        for feature in translation['protein_features']:
                            if feature.get('dbname', "") == "Pfam":
                                emitter.emit_edge(
                                    PFAMAlignment(start=int(feature['start']), end=int(feature['end'])),
                                    from_gid=prot.gid(),
                                    to_gid=PFAMFamily.make_gid(feature['name']))

emitter.close()
