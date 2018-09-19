
from bmeg.emitter import new_emitter
from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.vertex import Phenotype
from bmeg.edge import *  # noqa
from bmeg.enrichers.phenotype_enricher import normalize
from bmeg.gid import GID
from bmeg.stores import new_store

import glob
import logging
import sys
import ujson

DEFAULT_DIRECTORY = 'phenotype'


def transform(
    emitter_name="json",
    output_dir="outputs",
    emitter_directory=DEFAULT_DIRECTORY,
    vertex_names="**/*Phenotype.Vertex.json*",
    edge_names="**/*.Edge.json*",
    store_path="source/phenotype/sqlite.db"
):
    batch_size = 1000
    phenotype_cache = {}
    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix='normalized')
    path = '{}/{}'.format(output_dir, vertex_names)
    files = [filename for filename in glob.iglob(path, recursive=True) if 'normalized' not in filename]
    logging.info(files)
    logging.info(store_path)
    store = new_store('key-val', path=store_path, index=True)
    store.index()  # default is no index
    c = t = e = 0
    for file in files:
        logging.info(file)
        with reader(file) as ins:
            for line in ins:
                try:
                    # get the phenotype the transformer wrote
                    phenotype = ujson.loads(line)
                    phenotype_gid = phenotype['gid']
                    phenotype = phenotype['data']
                    # if un-normalized, normalize it
                    if phenotype['term'] == 'TODO':
                        # do we have it already?
                        stored_phenotype = store.get(phenotype['name'])
                        if not stored_phenotype:
                            # nope, fetch it
                            ontology_terms = normalize(phenotype['name'])
                            if len(ontology_terms) == 0:
                                # no hits? set term and id to name
                                phenotype['term'] = phenotype['name']
                                phenotype['term_id'] = 'NO_ONTOLOGY~{}'.format(phenotype['term'])
                            else:
                                # hits: set term and id to normalized term
                                phenotype['term'] = ontology_terms[0]['label']
                                phenotype['term_id'] = ontology_terms[0]['ontology_term']
                            # save it for next time
                            store.put(phenotype['name'], phenotype)
                        else:
                            phenotype = stored_phenotype
                    else:
                        # we have a phenotype with a term already
                        store.put(phenotype.get('name', phenotype.get('term')), phenotype)
                    phenotype = Phenotype.from_dict(phenotype)
                    emitter.emit_vertex(phenotype)
                    phenotype_cache[phenotype_gid] = phenotype
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
    files = [filename for filename in glob.iglob(path, recursive=True) if 'normalized' not in filename]
    c = t = e = 0
    for file in files:
        logging.info(file)
        with reader(file) as ins:
            for line in ins:
                try:
                    edge = ujson.loads(line)
                    if 'Phenotype:' not in edge['gid']:
                        logging.info('Edge {} has no phenotypes that need transformation. skipping.'.format(file))
                        break
                    # get edge components
                    label = edge['label']
                    from_ = edge['from']
                    to = edge['to']
                    data = edge['data']
                    # replace with normalized phenotype
                    if to in phenotype_cache:
                        to = phenotype_cache[to].gid()
                    if from_ in phenotype_cache:
                        from_ = phenotype_cache[from_].gid()
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
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)

    transform()
