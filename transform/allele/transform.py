
"""
dedupe and Allele.Vertex.json file into
"""
# import bmeg.enrichers.gene_enricher as gene_enricher
import logging
import glob
import sys
import json
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging

from bmeg.vertex import Allele
# from bmeg.edge import AlleleIn

""" """
store = {}


def get_allele(gid):
    return store.get(gid, None)


def put_allele(allele):
    store[allele.gid()] = allele


def dedupe(allele):
    existing_allele = get_allele(allele.gid())
    logging.debug('existing_allele {}'.format(existing_allele))
    return allele


def transform(output_dir, prefix, emitter_name='json'):
    path = '{}/**/*.Allele.Vertex.json'.format(output_dir)
    for filename in glob.iglob(path, recursive=True):
        logging.debug(filename)
        with open(filename, "r") as ins:
            for line in ins:
                allele = Allele.from_dict(json.loads(line)['data'])
                put_allele(dedupe(allele))


def main():  # pragma: no cover
    parser = default_argument_parser()
    parser.add_argument('--output_dir', type=str,
                        help='Path to the directory containing **/*.Allele.Vertex.json',
                        default='outputs')
    # We don't need the first argument, which is the program name
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform(mafpath=options.output_dir,
              prefix=options.prefix,
              emitter_name=options.emitter)


if __name__ == '__main__':  # pragma: no cover
    main()
