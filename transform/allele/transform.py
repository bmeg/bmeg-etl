
"""
dedupe and Allele.Vertex.json file into store
"""
import logging
import glob
import sys
import ujson
import dataclasses
import subprocess
import os.path
import threading

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


def enrich_from_dict(allele, myvariantinfo_annotation):
    """ assign or lookup """
    if not myvariantinfo_annotation:
        # no, fetch it from biothings
        myvariantinfo_annotation = myvariantinfo(
            genome=allele.genome, chromosome=allele.chromosome, start=allele.start,
            end=allele.end, reference_bases=allele.reference_bases,
            alternate_bases=allele.alternate_bases
        )
        if myvariantinfo_annotation:
            logging.debug('found new myvariant annotation for {}'.format(allele.gid()))
    allele.annotations.myvariantinfo = myvariantinfo_annotation
    return allele


def from_dict(allele_dicts):
    """ turn array of dicts into array of Alleles"""
    return [Allele.from_dict(d['data']) for d in allele_dicts]


def group_sorted_alleles(sorted_allele_file):
    """ yield an array of data records with the same gid"""
    with reader(sorted_allele_file) as ins:
        data = ujson.loads(ins.readline())
        _id = data['_id']
        alleles = [data]
        for line in ins:
            data = ujson.loads(line)
            if data['_id'] != _id:
                yield from_dict(alleles)
                alleles = []
            _id = data['_id']
            alleles.append(data)
        yield from_dict(alleles)


def sort_allele_files(path, sorted_allele_file):
    """ sort alleles file_names[] into a file in sorted_allele_file"""
    try:
        logging.debug(path)
        if not os.path.isfile(sorted_allele_file):
            files = [filename for filename in glob.iglob(path, recursive=True)]
            assert len(files) > 0
            files = ' '.join(files)
            cat = 'cat'
            if '.gz' in files:
                cat = 'zcat'
            # cat the files, sort on first key and save
            cmd = '{} {} | sort -k1 > {}'.format(cat, files, sorted_allele_file)
            logging.info('running {}'.format(cmd))
            subprocess.check_output(cmd, shell=True)
        return sorted_allele_file
    except subprocess.CalledProcessError as sort_error:
        raise ValueError("sort error code {} {}".format(sort_error.returncode, sort_error.output))


def transform(output_dir,
              prefix,
              allele_store_name='memory',
              allele_store_path='/tmp/sqlite.db',
              emitter_name='json',
              vertex_filename_pattern='**/*.Allele.Vertex.json.gz',
              myvariantinfo_path=None,
              sorted_allele_file='/tmp/sorted_allele_file.json'):
    """
    dedupe & reconcile annotations.
    * if sorted_allele_file does not exist create it by sorting all AlleleFiles decending from output_dir
    * reconcile AlleleAnnotations
    * enrich w/ myvariantinfo
    """
    threading.local().skip_check_types = True
    emitter = new_emitter(name=emitter_name, prefix=prefix)

    path = '{}/{}'.format(output_dir, vertex_filename_pattern)
    logging.info('checking {}'.format(path))
    sorted_allele_file = sort_allele_files(path, sorted_allele_file)
    allele_store = new_store('allele-sqlite', path=allele_store_path)
    c = 0
    t = 0
    batch_size = 20000
    batch = []
    gid_cache = {}
    if allele_store.size() == 0:
        logging.info('merging')
        for alleles in group_sorted_alleles(sorted_allele_file):
            allele = merge(alleles)
            gid = allele.gid()
            batch.append(
                (
                    gid,
                    allele.__class__.__name__,
                    ujson.dumps(dataclasses.asdict(allele))
                )
            )
            gid_cache[gid] = True
            c += 1
            t += 1
            if c % batch_size == 0:
                logging.info(t)
                allele_store.load_many(batch)
                c = 0
                batch = []
        logging.info(t)
        allele_store.load_many(batch)
        allele_store.index()
    else:
        logging.info('already merged, loading gid_cache')
        c = t = 0
        for gid in allele_store.all_ids():
            gid_cache[gid] = True
            c += 1
            t += 1
            if c % batch_size == 0:
                logging.info('loading gid_cache {}'.format(t))
                c = 0
        logging.info('loaded gid_cache {}'.format(t))

    logging.info('enriching')
    with reader(myvariantinfo_path) as ins:
        c = t = w = 0
        for line in ins:
            c += 1
            t += 1
            if c % batch_size == 0:
                logging.info('enriching read: {}, written: {}'.format(t, w))
                c = 0
            myvariant = {}
            try:
                myvariant = ujson.loads(line)
            except Exception as e:
                logging.warning(str(e))
            if 'hg19' not in myvariant:
                continue
            if 'vcf' not in myvariant:
                continue
            gid = Allele.make_gid(
                genome='GRCh37',
                chromosome=myvariant['chrom'],
                start=myvariant['hg19']['start'],
                end=myvariant['hg19']['end'],
                reference_bases=myvariant['vcf']['ref'],
                alternate_bases=myvariant['vcf']['alt'],
            )
            gid = str(gid)
            if gid not in gid_cache:
                continue
            allele = allele_store.get(gid)
            allele = enrich_from_dict(allele, myvariant)
            # TODO - is it worth batching these writes?
            allele_store.put(allele)
            w += 1
    logging.info('enriching finished read: {}, written: {}'.format(t, w))
    logging.info('enriching finished allele_store.size: {}'.format(allele_store.size()))
    logging.info('emitting')
    c = t = m = 0
    for allele in allele_store.all():
        if not allele.annotations.myvariantinfo:
            m += 1
        emitter.emit_vertex(allele)
        c += 1
        t += 1
        if c % batch_size == 0:
            logging.info('emitting written: {} myvariant misses: {}'.format(t, m))
            c = 0

    logging.info('emitting written: {} myvariant misses: {}'.format(t, m))
    emitter.close()


def main():  # pragma: no cover
    parser = default_argument_parser(prefix_default='allele')

    # TODO these should be exclusive
    parser.add_argument('--output_dir', type=str,
                        help='Path to the directory containing **/*.Allele.Vertex.json',
                        default='outputs')
    parser.add_argument('--sorted_allele_file', type=str,
                        help='Path to the file containing sorted **/*.Allele.Vertex.json',
                        default='source/allele/sorted_allele_file.json')
    parser.add_argument('--allele_store_name', type=str,
                        help='name of allele database [memory, sqlite, (mongo, bmeg)]',
                        default='allele-sqlite')
    parser.add_argument('--allele_store_path', type=str,
                        help='path of allele store',
                        default='source/allele/sqlite.db')

    parser.add_argument('--myvariantinfo_path', type=str,
                        help='path to myvariantinfo json',
                        default='source/myvariant.info/biothings_current_old_hg19.json.gz')

    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform(output_dir=options.output_dir,
              prefix=options.prefix,
              sorted_allele_file=options.sorted_allele_file,
              allele_store_name=options.allele_store_name,
              allele_store_path=options.allele_store_path,
              myvariantinfo_path=options.myvariantinfo_path,
              emitter_name=options.emitter)


if __name__ == '__main__':  # pragma: no cover
    main()
