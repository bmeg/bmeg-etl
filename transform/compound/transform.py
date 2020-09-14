#!/usr/bin/env python

import json
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


class MappingTables:
    def __init__(self, table_dir):
        self.table_dir = table_dir
        self.tables = {}

    def get(self, id_source, name):
        if id_source not in self.tables:
            d = {}
            tpath = os.path.join(self.table_dir, id_source + ".table")
            print("Opening Translation Table %s" % (tpath))
            with open(tpath) as handle:
                for line in handle:
                    row = line.rstrip().split("\t")
                    d[row[0]] = row[1]
            self.tables[id_source] = d
        return self.tables[id_source].get(name, None)

def transform(vertex_names="**/*Compound.Vertex.json*",
              edge_names="**/*Compound*.Edge.json*",
              output_dir="outputs",
              emitter_name="json",
              emitter_directory="compound",
              table_dir="reference/compound",
              biothings_source="source/compound/biothings.json"):

    missing_dir = "./missing"

    missing_handle = open(os.path.join(missing_dir, "compounds.tsv"), "w")

    mt = MappingTables(table_dir)

    biothings = {}
    with open(biothings_source) as handle:
        for line in handle:
            data = json.loads(line)
            if data is not None:
                id = data["id"]
                data['id'] = Compound.make_gid(id)
                compound = Compound(
                    project_id=Project.make_gid('Reference'),
                    **data
                )
                biothings[id] = compound


    batch_size = 1000
    compound_cache = {}
    dups = {}
    emitter = new_emitter(name=emitter_name, directory=emitter_directory, prefix='normalized')
    path = '{}/{}'.format(output_dir, vertex_names)
    vertex_files = [filename for filename in glob.iglob(path, recursive=True) if 'normalized' not in filename]
    logging.info(vertex_files)
    c = t = e = 0

    for file in vertex_files:
        logging.info(file)
        print("Parsing", file)
        with reader(file) as ins:
            for line in ins:
                # get the compound the transformer wrote
                compound = ujson.loads(line)
                compound_gid = compound['gid']
                compound_name = compound['data']['submitter_id']
                id_source = compound['data']['id_source']

                id = mt.get(id_source, compound_name)
                if id is None:
                    print("Compound %s:%s not found" % (id_source, compound_name))
                    missing_handle.write("%s\t%s\n" % (id_source, compound_name))
                    c = Compound(
                        id = Compound.make_gid("NOTFOUND:%s" % (compound_name)),
                        project_id=Project.make_gid('Reference'),
                        id_source=id_source
                    )
                    emitter.emit_vertex(c)
                    compound_cache[compound_gid] = c
                else:
                    print("Translate %s to %s" % (compound_gid, id))
                    emitter.emit_vertex(biothings[id])
                    compound_cache[compound_gid] = biothings[id]

    missing_handle.close()

    # get the edges
    path = '{}/{}'.format(output_dir, edge_names)
    edge_files = [filename for filename in glob.iglob(path, recursive=True) if 'normalized' not in filename]
    c = t = e = 0
    for file in edge_files:
        logging.info(file)
        print("Parsing %s" % (file))
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
