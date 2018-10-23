
from bmeg.emitter import new_emitter
from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.vertex import Publication
from bmeg.edge import *  # noqa

import glob
import logging
import sys
import ujson

DEFAULT_DIRECTORY = 'publication'


def transform(
    emitter_name="json",
    output_dir="outputs",
    emitter_directory=DEFAULT_DIRECTORY,
    edge_names="**/*.Edge.json*"
):
    batch_size = 1000
    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix='stub')

    # get the edges
    path = '{}/{}'.format(output_dir, edge_names)
    files = [filename for filename in glob.iglob(path, recursive=True) if 'stub' not in filename]
    c = t = e = 0
    for file in files:
        logging.info(file)
        with reader(file) as ins:
            for line in ins:
                try:
                    edge = ujson.loads(line)
                    if 'Publication:' not in edge['gid']:
                        logging.info('Edge {} has no publications that need transformation. skipping.'.format(file))
                        break
                    # get edge components
                    to = edge['to']
                    id = to.replace('Publication:', '')
                    title = abstract = text = date = author = citation = None
                    publication = Publication(id, title, abstract, text, date, author, citation)
                    emitter.emit_vertex(publication)
                    c += 1
                    t += 1
                except Exception as exc:
                    logging.warning(str(exc))
                    raise
                    e += 1
                if c % batch_size == 0:
                    logging.info('transforming read: {} errors: {}'.format(t, e))
                    c = 0
        logging.info('transforming read: {} errors: {}'.format(t, e))
        emitter.close()


if __name__ == '__main__':  # pragma: no cover
    parser = default_argument_parser()
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform()
