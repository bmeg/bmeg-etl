from bmeg.emitter import new_emitter
from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.enrichers.phenotype_enricher import normalize
from bmeg.stores import new_store
from bmeg import *  # noqa: F403

import logging
import glob
import os
import re
import sys
import ujson


def transform(
    vertex_names="**/*Phenotype.Vertex.json*",
    edge_names="**/*Phenotype*.Edge.json*",
    output_dir="outputs",
    emitter_name="json",
    emitter_directory="phenotype",
    store_path="source/phenotype/sqlite.db"
):
    batch_size = 1000
    phenotype_cache = {}
    dups = {}
    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix='normalized')

    path = '{}/{}'.format(output_dir, vertex_names)
    vertex_files = [filename for filename in glob.iglob(path, recursive=True) if 'normalized' not in filename and 'mondo' not in filename]
    logging.info("vertex files: %s", vertex_files)

    logging.info("store path: %s", store_path)
    store = new_store('key-val', path=store_path, index=True)
    store.index()  # default is no index

    c = t = e = 0
    for file in vertex_files:
        logging.info("processing file: %s", file)
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
                            ontology_term = normalize(phenotype['name'])
                            if ontology_term is None:
                                # no hits? set term and id to name
                                phenotype['term'] = phenotype['name']
                                phenotype['term_id'] = 'NO_ONTOLOGY:{}'.format(phenotype['term'])
                                phenotype['id'] = Phenotype.make_gid('NO_ONTOLOGY:{}'.format(phenotype['term']))
                            else:
                                # hits: set term and id to normalized term
                                phenotype['term'] = ontology_term['label']
                                phenotype['term_id'] = ontology_term['ontology_term']
                                phenotype['id'] = Phenotype.make_gid(ontology_term['ontology_term'])
                            # save it for next time
                            store.put(phenotype['name'], phenotype)
                        else:
                            phenotype = stored_phenotype
                            phenotype['id'] = Phenotype.make_gid(phenotype['term_id'])
                    elif not phenotype['term_id'].startswith('MONDO'):
                        stored_phenotype = store.get(phenotype['term_id'])
                        if not stored_phenotype:
                            # nope, fetch it
                            ontology_term = normalize(phenotype['term_id'])
                            if ontology_term is None:
                                ontology_term = normalize(phenotype['term'])
                            if ontology_term is None:
                                phenotype['term'] = phenotype['term']
                                phenotype['term_id'] = 'NO_ONTOLOGY:{}'.format(phenotype['term_id'])
                                phenotype['id'] = Phenotype.make_gid('NO_ONTOLOGY:{}'.format(phenotype['term_id']))
                            else:
                                # hits: set term and id to normalized term
                                phenotype['term'] = ontology_term['label']
                                phenotype['term_id'] = ontology_term['ontology_term']
                                phenotype['id'] = Phenotype.make_gid(ontology_term['ontology_term'])
                            # save it for next time
                            store.put(phenotype['term_id'], phenotype)
                        else:
                            phenotype = stored_phenotype
                            phenotype['id'] = Phenotype.make_gid(phenotype['term_id'])
                    else:
                        # we have a phenotype with a term already
                        phenotype['id'] = Phenotype.make_gid(phenotype['term_id'])
                        store.put(phenotype['term_id'], phenotype)
                        store.put(phenotype['term'], phenotype)

                    phenotype = Phenotype(**phenotype)
                    if phenotype.gid() not in dups:
                        emitter.emit_vertex(phenotype)
                        dups[phenotype.gid()] = None
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
    edge_files = [filename for filename in glob.iglob(path, recursive=True) if 'normalized' not in filename and 'mondo' not in filename]
    logging.info("edge files: %s", edge_files)

    c = t = e = 0
    for file in edge_files:
        logging.info("processing file: %s", file)
        with reader(file) as ins:
            for line in ins:
                try:
                    edge = ujson.loads(line)
                    if 'Phenotype:' not in edge['gid']:
                        logging.info('Edge {} has no phenotypes that need transformation. skipping.'.format(file))
                        break

                    # get edge components
                    from_ = edge['from']
                    to = edge['to']
                    data = edge['data']

                    # replace with normalized phenotype
                    if to in phenotype_cache:
                        to = phenotype_cache[to].gid()
                    if from_ in phenotype_cache:
                        from_ = phenotype_cache[from_].gid()

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
