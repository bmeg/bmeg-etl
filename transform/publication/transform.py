
from bmeg.emitter import new_emitter
from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg import Publication, Project

import glob
import logging
import sys
import ujson


def transform(
        emitter_name="json",
        output_dir="outputs",
        emitter_directory="publication",
        vertex_names="**/*Publication.Vertex.json*",
        edge_names="**/*Publication*.Edge.json*"
):
    # get the existing publication vertices
    path = '{}/{}'.format(output_dir, vertex_names)
    files = [filename for filename in glob.iglob(path, recursive=True) if 'stub' not in filename]
    nfiles = len(files)
    dedup = {}
    f = v = 0
    for file in files:
        f += 1
        logging.info("processing publication file: {}/{}".format(f, nfiles))
        with reader(file) as ins:
            for line in ins:
                try:
                    vertex = ujson.loads(line)
                    if 'Publication' != vertex['label']:
                        logging.info('Vertex {} does not appear to contain publications. skipping.'.format(file))
                        break
                    dedup[vertex["gid"]] = True
                    v += 1
                except Exception as exc:
                    logging.error(str(exc))
                    raise exc
    logging.info('processed publication vertices: {}'.format(v))

    # get the edges connected to publications
    batch_size = 1000
    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix='stub')
    path = '{}/{}'.format(output_dir, edge_names)
    files = [filename for filename in glob.iglob(path, recursive=True) if 'stub' not in filename]
    nfiles = len(files)
    e = f = r = 0
    for file in files:
        f += 1
        logging.info("processing HasSupportingReference file: {}/{}".format(f, nfiles))
        with reader(file) as ins:
            for line in ins:
                try:
                    edge = ujson.loads(line)
                    if 'Publication:' not in edge['gid']:
                        logging.info('Edge {} has no publications that need transformation. skipping.'.format(file))
                        break
                    # get edge components
                    to = edge['to']
                    if to in dedup:
                        r += 1
                        continue
                    dedup[to] = True
                    url = to.replace('Publication:', 'http://')
                    publication = Publication(
                        id=Publication.make_gid(url),
                        url=url,
                        project_id=Project.make_gid("Reference")
                    )
                    emitter.emit_vertex(publication)
                    e += 1
                except Exception as exc:
                    logging.error(str(exc))
                    raise exc
                if e % batch_size == 0:
                    logging.info('emitted stub publication vertices: {}'.format(e))
    logging.info('emitted stub publication vertices: {}'.format(e))
    logging.info('existing publication refs found: {}'.format(r))
    emitter.close()


if __name__ == '__main__':  # pragma: no cover
    parser = default_argument_parser()
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)
    transform()
