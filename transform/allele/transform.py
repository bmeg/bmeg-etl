
"""
dedupe and Allele.Vertex.json file into store
"""
import logging
import glob
import sys
import json
import uuid
import dataclasses
import subprocess
import os.path

from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.emitter import new_emitter
from bmeg.vertex import Allele, AlleleAnnotations
from bmeg.enrichers.allele_enricher import myvariantinfo
from bmeg.stores import new_store
from bmeg.ioutils import reader


def merge(alleles):
    """ reconcile annotations """
    # first one in list is reference
    alleles = iter(alleles)
    base = next(alleles)
    for c in alleles:
        for f in dataclasses.fields(AlleleAnnotations):
            existing = getattr(base.annotations, f.name)
            current = getattr(c.annotations, f.name)
            if not existing and current:
                setattr(base.annotations, f.name, current)
    return base


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


def from_dict(allele_dicts):
    """ turn array of dicts into array of Alleles"""
    return [Allele.from_dict(d['data']) for d in allele_dicts]

def group_sorted_alleles(sorted_allele_file):
    """ yield an array of data records with the same gid"""
    with reader(sorted_allele_file) as ins:
        data = json.loads(ins.readline())
        _id = data['_id']
        alleles = [data]
        for line in ins:
            data = json.loads(line)
            if data['_id'] != _id:
                yield from_dict(alleles)
                alleles = []
            _id = data['_id']
            alleles.append(data)
        yield from_dict(alleles)


def sort_allele_files(path, sorted_alleles_file):
    """ sort alleles file_names[] into a file in sorted_alleles_file"""
    try:
        logging.debug(path)
        if not os.path.isfile(sorted_alleles_file):
            files = [filename for filename in glob.iglob(path, recursive=True)]
            assert len(files) > 0
            files = ' '.join(files)
            cat = 'cat'
            if '.gz' in files:
                cat = 'zcat'
            # cat the files, sort on first key and save
            cmd = '{} {} | sort -k1 > {}'.format(cat, files, sorted_alleles_file)
            logging.info('running {}'.format(cmd))
            subprocess.check_output(cmd, shell=True)
        return sorted_alleles_file
    except subprocess.CalledProcessError as sort_error:
        raise ValueError("sort error code {} {}".format(sort_error.returncode, sort_error.output))


def transform(output_dir, prefix,
              myvariant_store_name='memory',
              emitter_name='json',
              vertex_filename_pattern='**/*.Allele.Vertex.json.gz',
              myvariantinfo_path=None,
              sorted_alleles_file='/tmp/sorted_allele_file.json'):
    """
    dedupe & reconcile annotations.
    * if sorted_allele_file does not exist create it by sorting all AlleleFiles decending from output_dir
    * reconcile AlleleAnnotations
    * enrich w/ myvariantinfo
    """

    emitter = new_emitter(name=emitter_name, prefix=prefix)

    logging.info('creating myvariant_store')
    myvariant_store = new_store(myvariant_store_name, myvariantinfo_path=myvariantinfo_path)
    c = 0
    t = 0
    batch_size = 20000
    logging.info('sorting')
    path = '{}/{}'.format(output_dir, vertex_filename_pattern)
    sorted_allele_file = sort_allele_files(path, sorted_alleles_file)
    logging.info('merging/enrich')
    for alleles in group_sorted_alleles(sorted_allele_file):
        allele = merge(alleles)
        allele = enrich(allele, myvariant_store)
        emitter.emit_vertex(allele)
        c += 1
        t += 1
        if c % batch_size == 0:
            logging.info(t)
            c = 0
    emitter.close()


def main():  # pragma: no cover
    parser = default_argument_parser()

    # TODO these should be exclusive
    parser.add_argument('--output_dir', type=str,
                        help='Path to the directory containing **/*.Allele.Vertex.json',
                        default='outputs')
    parser.add_argument('--sorted_allele_file', type=str,
                        help='Path to the file containing sorted **/*.Allele.Vertex.json',
                        default='/tmp/sorted_allele_file.json')

    parser.add_argument('--myvariant_store_name', type=str,
                        help='name of allele database [memory, sqlite, mongo, bmeg]',
                        default='myvariantinfo-memory')
    parser.add_argument('--myvariantinfo_path', type=str,
                        help='path to myvariantinfo json',
                        default='source/myvariant.info/biothings_current_old_hg19.json.gz')

    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform(output_dir=options.output_dir,
              prefix=options.prefix,
              myvariant_store_name=options.myvariant_store_name,
              myvariantinfo_path=options.myvariantinfo_path,
              emitter_name=options.emitter)


if __name__ == '__main__':  # pragma: no cover
    main()
