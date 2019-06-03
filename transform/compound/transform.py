
from bmeg.emitter import new_emitter
from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.vertex import Compound
from bmeg.edge import *  # noqa
from bmeg.enrichers.drug_enricher import normalize
from bmeg.gid import GID
from bmeg.stores import new_store

import logging
import sys
import ujson

DEFAULT_DIRECTORY = 'compound'


def transform(vertex_files, edge_files,
              emitter_name="json",
              emitter_directory=DEFAULT_DIRECTORY,
              store_path="source/compound/sqlite.db"):
    batch_size = 1000
    compound_cache = {}
    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix='normalized')
    logging.info(vertex_files)
    store = new_store('key-val', path=store_path, index=True)
    c = t = e = 0
    dups = []
    for file in vertex_files:
        logging.info(file)
        with reader(file) as ins:
            for line in ins:
                try:
                    # get the compound the transformer wrote
                    compound = ujson.loads(line)
                    compound_gid = compound['gid']
                    compound = compound['data']
                    # if un-normalized, normalize it
                    if compound['term'] == 'TODO':
                        # do we have it already?
                        stored_compound = store.get(compound['name'])
                        if not stored_compound:
                            # nope, fetch it
                            ontology_terms = normalize(compound['name'])
                            if len(ontology_terms) == 0:
                                # no hits? set term and id to name
                                compound['term'] = compound['name']
                                compound['term_id'] = 'NO_ONTOLOGY~{}'.format(compound['term'])
                            else:
                                # hits: set term and id to normalized term
                                compound['term'] = ontology_terms[0]['synonym']
                                compound['term_id'] = ontology_terms[0]['ontology_term']
                            # save it for next time
                            store.put(compound['name'], compound)
                        else:
                            compound = stored_compound
                    else:
                        # we have a compound with a term already
                        store.put(compound['name'], compound)
                    compound = Compound.from_dict(compound)
                    if compound.gid() not in dups:
                        emitter.emit_vertex(compound)
                    compound_cache[compound_gid] = compound
                    dups.append(compound.gid())
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

    # get the edges
    c = t = e = 0
    for file in edge_files:
        logging.info(file)
        with reader(file) as ins:
            for line in ins:
                try:
                    edge = ujson.loads(line)
                    if 'Compound:' not in edge['gid']:
                        logging.info('Edge {} has no compounds that need transformation. skipping.'.format(file))
                        break
                    # get edge components
                    label = edge['label']
                    from_ = edge['from']
                    to = edge['to']
                    data = edge['data']
                    # replace with normalized compound
                    if to in compound_cache:
                        to = compound_cache[to].gid()
                    if from_ in compound_cache:
                        from_ = compound_cache[from_].gid()
                    cls = globals()[label]
                    edge = cls()
                    if data:
                        edge = cls(data=data)
                    emitter.emit_edge(
                        edge,
                        GID(from_),
                        GID(to),
                    )
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
    parser.add_argument('--vertex_files', nargs='+', help='vertex list', required=True)
    parser.add_argument('--edge_files', nargs='+', help='edge list', required=True)
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform(vertex_files=options.vertex_files, edge_files=options.edge_files)
