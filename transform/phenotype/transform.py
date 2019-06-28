from bmeg import Phenotype, GenericEdge
from bmeg.emitter import new_emitter
from bmeg.ioutils import reader
from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging
from bmeg.enrichers.phenotype_enricher import normalize
from bmeg.stores import new_store

import logging
import sys
import ujson


def transform(
    vertex_files,
    edge_files,
    emitter_name="json",
    emitter_directory="phenotype",
    store_path="source/phenotype/sqlite.db"
):
    batch_size = 1000
    phenotype_cache = {}
    dups = {}
    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix='normalized')

    logging.info("vertex files:", edge_files)
    logging.info(store_path)
    store = new_store('key-val', path=store_path, index=True)
    store.index()  # default is no index
    c = t = e = 0
    for file in vertex_files:
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
                                phenotype['id'] = Phenotype.make_gid('NO_ONTOLOGY~{}'.format(phenotype['term']))
                            else:
                                # hits: set term and id to normalized term
                                phenotype['term'] = ontology_terms[0]['label']
                                phenotype['term_id'] = ontology_terms[0]['ontology_term']
                                phenotype['id'] = Phenotype.make_gid(ontology_terms[0]['ontology_term'])
                            # save it for next time
                            store.put(phenotype['name'], phenotype)
                        else:
                            phenotype = stored_phenotype
                            phenotype['id'] = Phenotype.make_gid(phenotype['term_id'])
                    else:
                        if 'MONDO' not in phenotype['term_id']:
                            # we prefer MONDO
                            # do we have it already?
                            stored_phenotype = store.get(phenotype['term'])
                            if not stored_phenotype:
                                # nope, fetch it
                                ontology_terms = normalize(phenotype['term'])
                                if len(ontology_terms) > 0:
                                    # hits: set term and id to normalized term
                                    phenotype['term'] = ontology_terms[0]['label']
                                    phenotype['term_id'] = ontology_terms[0]['ontology_term']
                                    phenotype['id'] = Phenotype.make_gid(ontology_terms[0]['ontology_term'])

                        # we have a phenotype with a term already
                        phenotype['id'] = Phenotype.make_gid(phenotype['term_id'])
                        store.put(phenotype.get('name', phenotype.get('term')), phenotype)

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
    logging.info("edge files:", edge_files)
    c = t = e = 0
    for file in edge_files:
        logging.info(file)
        with reader(file) as ins:
            for line in ins:
                try:
                    edge = ujson.loads(line)
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

                    emitter.emit_edge(
                        GenericEdge(
                            from_gid=from_,
                            to_gid=to,
                            edge_label=label,
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
    parser.add_argument('--vertex_files', nargs='+', help='vertex list', required=True)
    parser.add_argument('--edge_files', nargs='+', help='edge list', required=True)
    options = parser.parse_args(sys.argv[1:])
    default_logging(options.loglevel)
    transform(vertex_files=options.vertex_files, edge_files=options.edge_files)
