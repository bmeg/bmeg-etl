
"""
dedupe and Allele.Vertex.json file into store
"""
import logging
import glob
import sys
import json
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.emitter import new_emitter
from bmeg.vertex import Allele, AlleleAnnotations
from bmeg.enrichers.allele_enricher import myvariantinfo
from bmeg.stores import new_store

import dataclasses


def merge(allele, existing_allele):
    """ reconcile annotations """
    if existing_allele:
        for f in dataclasses.fields(AlleleAnnotations):
            existing = getattr(existing_allele.annotations, f.name)
            current = getattr(allele.annotations, f.name)
            if not current and existing:
                setattr(allele, f.name, existing)
    return allele


def enrich(allele, myvariant_store):
    """ check store, fetch & update.  TODO - other enrichments?"""
    # skip if already done
    if not allele.annotations.myvariantinfo:
        # in local store?
        myvariantinfo_annotation = myvariant_store.get(allele.gid())
        if not myvariantinfo_annotation:
            # no, fetch it from biothings
            myvariantinfo_annotation = myvariantinfo(
                genome=allele.genome, chromosome=allele.chromosome, start=allele.start,
                end=allele.end, reference_bases=allele.reference_bases,
                alternate_bases=allele.alternate_bases
            )
            if myvariantinfo_annotation:
                logging.debug('found new myvariant annotation for {}'.format(allele.gid()))
                # save it for next time
                myvariant_store.put(allele)
        allele.annotations.myvariantinfo = myvariantinfo_annotation
    return allele


def transform(output_dir, prefix,
              store_name='memory',
              myvariant_store_name='memory',
              emitter_name='json',
              vertex_filename_pattern='**/*.Allele.Vertex.json',
              myvariantinfo_path=None):
    """
    dedupe & reconcile annotations.
    * look at all AllelFiles decending from output_dir
    * reconcile AlleleAnnotations
    * enrich w/ myvariantinfo
    """
    # TODO pass pattern, don't hardcode
    path = '{}/{}'.format(output_dir, vertex_filename_pattern)
    store = new_store(store_name)
    emitter = new_emitter(name=emitter_name, prefix=prefix)

    for filename in glob.iglob(path, recursive=True):
        with open(filename, "r") as ins:
            for line in ins:
                allele = Allele.from_dict(json.loads(line)['data'])
                # allele's nulls set tot store's alleles values
                store.put(merge(allele, store.get(allele.gid())))

    myvariant_store = new_store(myvariant_store_name, myvariantinfo_path=myvariantinfo_path)
    for allele in store.all():
        allele = enrich(allele, myvariant_store)
        emitter.emit_vertex(allele)

    emitter.close()


def main():  # pragma: no cover
    parser = default_argument_parser()
    # TODO should this just be wildcard
    parser.add_argument('--output_dir', type=str,
                        help='Path to the directory containing **/*.Allele.Vertex.json',
                        default='outputs')
    parser.add_argument('--allele_store_name', type=str,
                        help='name of allele database [memory, sqlite, mongo, bmeg]',
                        default='memory')
    parser.add_argument('--myvariant_store_name', type=str,
                        help='name of allele database [memory, sqlite, mongo, bmeg]',
                        default='myvariantinfo-memory')

    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform(mafpath=options.output_dir,
              prefix=options.prefix,
              store_name=options.store_name,
              myvariant_store_name=options.myvariant_store_name,
              emitter_name=options.emitter)


if __name__ == '__main__':  # pragma: no cover
    main()
