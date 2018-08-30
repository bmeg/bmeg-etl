
"""
dedupe and Allele.Vertex.json file into STORE
"""
# import bmeg.enrichers.gene_enricher as gene_enricher
import logging
import glob
import sys
import json
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.emitter import new_emitter
from bmeg.vertex import Allele, AlleleAnnotations
import dataclasses

# publish the store so the unit test can read it
STORE = None

class MemoryStore:
    """ store and retrieve Alleles """
    def __init__(self):
        self.key_val = {}
    def get_allele(self, gid):
        return self.key_val.get(gid, None)

    def put_allele(self, allele):
        self.key_val[allele.gid()] = allele

    def get_all(self):
        for k in self.key_val.keys():
            yield self.key_val[k]

def new_store(name):
    if name == 'memory':
        return MemoryStore()
    assert False, 'no store named {}'.format(name)

def merge(allele, existing_allele):
    """ reconcile annotations """
    if existing_allele:
        for f in dataclasses.fields(AlleleAnnotations):
            existing = getattr(existing_allele.annotations, f.name)
            current = getattr(allele.annotations, f.name)
            if not current and existing:
                setattr(allele, f.name, existing)
    # logging.debug('existing_allele {}'.format(existing_allele))
    return allele


def transform(output_dir, prefix, store_name='memory', emitter_name='json'):
    """ dedupe reconcile annotations """
    path = '{}/**/*.Allele.Vertex.json'.format(output_dir)
    STORE = new_store(store_name)
    emitter = new_emitter(name=emitter_name, prefix=prefix)

    for filename in glob.iglob(path, recursive=True):
        with open(filename, "r") as ins:
            for line in ins:
                allele = Allele.from_dict(json.loads(line)['data'])
                STORE.put_allele(merge(allele, STORE.get_allele(allele.gid())))

    for allele in STORE.get_all():
        # TODO fetch myvariant info
        emitter.emit_vertex(allele)
    emitter.close()    



def main():  # pragma: no cover
    parser = default_argument_parser()
    parser.add_argument('--output_dir', type=str,
                        help='Path to the directory containing **/*.Allele.Vertex.json',
                        default='outputs')
    parser.add_argument('--store_name', type=str,
                        help='name of allele database [memory, sqlite, mongo, bmeg]',
                        default='memory')

    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform(mafpath=options.output_dir,
              prefix=options.prefix,
              store_name=options.store_name,
              emitter_name=options.emitter)


if __name__ == '__main__':  # pragma: no cover
    main()
