
from bmeg.emitter import new_emitter
from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.enrichers.drug_enricher import normalize
from bmeg.stores import new_store
from bmeg import *  # noqa: F403

import logging
import glob
import os
import re
import sys
import ujson


def transform(vertex_names="**/*Compound.Vertex.json*",
              edge_names="**/*Compound*.Edge.json*",
              output_dir="outputs",
              emitter_name="json",
              emitter_directory="compound",
              store_path="source/compound/sqlite.db"):

    batch_size = 1000
    compound_cache = {}
    dups = {}
    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix='normalized')
    path = '{}/{}'.format(output_dir, vertex_names)
    vertex_files = [filename for filename in glob.iglob(path, recursive=True) if 'normalized' not in filename]
    logging.info(vertex_files)
    store = new_store('key-val', path=store_path, index=True)
    c = t = e = 0
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
                    if compound['id_source'] == 'TODO':
                        # do we have it already?
                        stored_compound = store.get(compound['submitter_id'])
                        if not stored_compound:
                            # nope, fetch it
                            cinfo = normalize(compound['submitter_id'])
                            if cinfo is None:
                                # no hits? set term and id to name
                                compound['id_source'] = 'NO_ONTOLOGY'
                                compound['id'] = Compound.make_gid('NO_ONTOLOGY:{}'.format(compound['submitter_id']))
                            else:
                                compound.update(cinfo)
                                compound['id'] = Compound.make_gid(compound['id'])
                            # save it for next time
                            store.put(compound['submitter_id'], compound)
                        else:
                            compound = stored_compound
                            compound['id'] = Compound.make_gid(compound['id'].strip('Compound:'))
                    else:
                        # we have a compound with a term already
                        compound['id'] = Compound.make_gid(compound['id'].strip('Compound:'))
                        store.put(compound['submitter_id'], compound)

                    # create compound and emit
                    compound = Compound(**compound)
                    if compound.gid() not in dups:
                        emitter.emit_vertex(compound)
                        dups[compound.gid()] = None
                    compound_cache[compound_gid] = compound
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
    path = '{}/{}'.format(output_dir, edge_names)
    edge_files = [filename for filename in glob.iglob(path, recursive=True) if 'normalized' not in filename]
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
                    from_ = edge['from']
                    to = edge['to']
                    data = edge['data']

                    # replace with normalized compound if available
                    if to in compound_cache:
                        to = compound_cache[to].gid()
                    if from_ in compound_cache:
                        from_ = compound_cache[from_].gid()

                    etype = re.sub(".Edge.json.*", "", os.path.basename(file)).split(".")[-1]
                    from_label, edge_label, to_label = etype.split("_")
                    from_cls = globals()[from_label]
                    to_cls = globals()[to_label]
                    edge_cls = globals()[etype]

                    emitter.emit_edge(
                        edge_cls(
                            from_gid=from_cls._gid_cls(from_),
                            to_gid=to_cls._gid_cls(to),
                            data=data
                        )
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
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)
    transform()
