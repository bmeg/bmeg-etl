
"""
dedupe and Allele.Vertex.json file into store
"""
import logging
import glob
import sys
import ujson
import subprocess
import os.path
import threading

from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.emitter import new_emitter
from bmeg import Allele
from bmeg.ioutils import reader


def valid_filenames(path):
    """ filter out any we do not want to process """
    filenames = []
    for filename in glob.iglob(path, recursive=True):
        # myvariant download **DEPRECIATED**
        if 'myvariant' in filename:
            continue
        # our own output
        if 'allele.' in filename:
            continue
        # pseudo alleles **DEPRECIATED**
        if 'MinimalAllele' in filename:
            continue
        filenames.append(filename)
    return filenames


def sort_allele_files(path, sorted_allele_file):
    """ sort alleles file_names[] into a file in sorted_allele_file"""
    try:
        logging.info(path)
        logging.info(sorted_allele_file)
        if not os.path.isfile(sorted_allele_file):
            files = valid_filenames(path)
            assert len(files) > 0
            files = ' '.join(files)
            cat = 'cat'
            if '.gz' in files:
                cat = 'zcat'
            # cat the files, sort on first key and save
            cmd = '{} {} | sort -k1 | uniq  > {}'.format(cat, files, sorted_allele_file)
            logging.info('running {}'.format(cmd))
            subprocess.check_output(cmd, shell=True)
        return sorted_allele_file
    except subprocess.CalledProcessError as sort_error:
        raise ValueError("sort error code {} {}".format(sort_error.returncode, sort_error.output))


def transform(output_dir,
              prefix,
              emitter_name='json',
              vertex_filename_pattern='**/*Allele.Vertex.json.gz',
              myvariantinfo_path=None,
              sorted_allele_file='/tmp/sorted_allele_file.json'):
    """
    dedupe & reconcile annotations.
    * if sorted_allele_file does not exist create it by sorting all AlleleFiles decending from output_dir
    * reconcile AlleleAnnotations
    * enrich w/ myvariantinfo
    """
    threading.local().skip_check_types = True
    emitter = new_emitter(name=emitter_name, directory=prefix)
    # Step one:
    # * sort  and dedupe the output/**/*.Allele.Vertex.json.gz files
    path = '{}/{}'.format(output_dir, vertex_filename_pattern)
    logging.info('checking {}'.format(path))
    sorted_allele_file = sort_allele_files(path, sorted_allele_file)
    # * emit all alleles
    logging.info('emitting')
    c = t = 0
    batch_size = 20000
    with reader(sorted_allele_file) as ins:
        for line in ins:
            try:
                allele_dict = ujson.loads(line)['data']
                allele_dict['id'] = Allele.make_gid(
                    allele_dict['genome'], allele_dict['chromosome'],
                    allele_dict['start'], allele_dict['end'],
                    allele_dict['reference_bases'], allele_dict['alternate_bases']
                )
                allele = Allele(**allele_dict)
                emitter.emit_vertex(allele)
                c += 1
                t += 1
                if c % batch_size == 0:
                    logging.info('emitting written: {}'.format(t))
                    c = 0
            except Exception as e:
                print('{} {}'.format(line, str(e)))
    logging.info('emitting written: {}'.format(t))
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
    parser.add_argument('--delete_workfile', dest='delete_workfile', action='store_true')
    parser.add_argument('--no-delete_workfile', dest='delete_workfile', action='store_false')
    parser.set_defaults(delete_workfile=True)

    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    if options.delete_workfile and os.path.isfile(options.sorted_allele_file):
        logging.info('deleting {}'.format(options.sorted_allele_file))
        os.remove(options.sorted_allele_file)

    transform(output_dir=options.output_dir,
              prefix=options.prefix,
              sorted_allele_file=options.sorted_allele_file,
              emitter_name=options.emitter)


if __name__ == '__main__':  # pragma: no cover
    main()
